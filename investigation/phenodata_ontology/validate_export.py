"""Independently read exported TSVs back and validate them against source CSVs."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from . import catalog
from .pipeline import DEFAULT_DATA, DEFAULT_OUTPUT, digest, load_inputs


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check_terms(annotation, location, column=False):
    for key in ("term_id", "unit_id", "column_term_id"):
        identifier = annotation.get(key, "")
        if identifier:
            require(catalog.is_current(identifier), f"{location}: invalid/obsolete {identifier}")
            if key != "column_term_id" and not column:
                require(annotation.get("status") != "needs_review",
                        f"{location}: review candidate leaked into assigned {key}")
    assertions = annotation.get("assertions") or []
    if isinstance(assertions, str):
        assertions = json.loads(assertions)
    for assertion in assertions:
        check_terms(assertion, location + "/assertion")


def require_same_rows(expected, actual, location):
    """Compare parsed rows, independent of field order and TSV quoting."""
    fields = sorted({key for row in expected + actual for key in row})
    signature = lambda row: tuple(row.get(field, "") for field in fields)
    require(Counter(map(signature, expected)) == Counter(map(signature, actual)),
            f"Project export differs from aggregate: {location}")


def expected_role_coverage(cells, values):
    coverage = []
    for role in sorted({cell["semantic_role"] for cell in cells}):
        role_cells = [cell for cell in cells if cell["semantic_role"] == role]
        role_values = [value for value in values if value["semantic_role"] == role]
        coverage.append({"semantic_role": role, "cells": str(len(role_cells)),
                         "unique_value_mappings": str(len(role_values)),
                         "cells_with_assigned_term": str(sum(bool(cell.get("term_id")) for cell in role_cells)),
                         "values_with_assigned_term": str(sum(bool(value.get("term_id")) for value in role_values)),
                         "cells_with_assigned_unit": str(sum(bool(cell.get("unit_id")) for cell in role_cells)),
                         "values_with_assigned_unit": str(sum(bool(value.get("unit_id")) for value in role_values)),
                         "cells_needing_review": str(sum(cell["status"] == "needs_review" for cell in role_cells))})
    return coverage


def validate_project_exports(output, report, tables, project_output_root=None):
    """Check every per-project artifact against the aggregate source of truth."""
    expected_projects = {row["project"] for row in tables["inventory"]}
    exports = report.get("project_exports", {})
    if project_output_root is None and not exports:
        return 0
    if project_output_root is None:
        require(set(exports) == expected_projects, "Project export scope differs from aggregate inventory")
        destinations = {project: (output / exports[project]["directory"]).resolve()
                        for project in expected_projects}
    else:
        destinations = {project: Path(project_output_root) / project / "phenodata"
                        for project in expected_projects}
    aggregate_coverage = read_tsv(output / "coverage_by_project.tsv")
    for project, directory in sorted(destinations.items()):
        require(directory.is_dir(), f"Missing project export: {directory}")
        selected = {}
        for name, aggregate in tables.items():
            expected = [row for row in aggregate if row.get("project") == project]
            actual = read_tsv(directory / (name + ".tsv"))
            require(all(row.get("project") == project for row in actual),
                    f"Another project's rows leaked into {directory}/{name}.tsv")
            require_same_rows(expected, actual, f"{project}/{name}.tsv")
            selected[name] = expected
        require_same_rows([row for row in aggregate_coverage if row["project"] == project],
                          read_tsv(directory / "coverage_by_project.tsv"),
                          f"{project}/coverage_by_project.tsv")
        require_same_rows(expected_role_coverage(selected["cell_annotations"], selected["value_mappings"]),
                          read_tsv(directory / "coverage_by_role.tsv"),
                          f"{project}/coverage_by_role.tsv")
        provenance = json.loads((directory / "provenance.json").read_text())
        require(provenance.get("project") == project, f"Incorrect project provenance: {directory}")
        for key in ("rule_file_sha256", "catalog_file_sha256", "evidence_file_sha256"):
            require(provenance.get(key) == report.get(key), f"Project provenance hashes differ: {project}/{key}")
        expected_counts = {"files": len(selected["inventory"]),
                           "source_cells": len(selected["cell_annotations"]),
                           "source_records": len(selected["row_annotations"]),
                           "source_column_instances": len(selected["column_mappings"]),
                           "unique_value_mappings": len(selected["value_mappings"])}
        for key, count in expected_counts.items():
            require(provenance.get(key) == count, f"Project provenance count differs: {project}/{key}")
        reference = provenance.get("aggregate_report", "")
        require(bool(reference) and (directory / reference).resolve() == (output / "qa_report.json").resolve(),
                f"Broken aggregate report reference: {project}")
    return len(destinations)


def validate(output=DEFAULT_OUTPUT, data_root=DEFAULT_DATA, project_output_root=None):
    output, data_root = Path(output), Path(data_root)
    report = json.loads((output / "qa_report.json").read_text())
    scope = report["scope_projects"]
    inputs = load_inputs(data_root, set(scope) if isinstance(scope, list) else None)
    source = {}
    for item in inputs:
        for record in item["records"]:
            for index, value in enumerate(record["cells"], 1):
                source[(item["relative_path"], record["source_record"], index)] = {
                    "value": value, "header": item["headers"][index - 1],
                    "sha256": item["sha256"], "line_start": record["source_line_start"],
                    "line_end": record["source_line_end"]}
    cells = read_tsv(output / "cell_annotations.tsv")
    values = read_tsv(output / "value_mappings.tsv")
    columns = read_tsv(output / "column_mappings.tsv")
    rows = read_tsv(output / "row_annotations.tsv")
    inventory = read_tsv(output / "inventory.tsv")
    curation = read_tsv(output / "curation_queue.tsv")
    require(len(cells) == len(source) == report["source_cells"], "Source cell coverage differs after TSV serialization")
    require(len(rows) == report["source_records"], "Row annotation count mismatch")
    require(len(columns) == report["source_column_instances"], "Column annotation count mismatch")
    require(len(inventory) == len(inputs) == report["files"], "Inventory count mismatch")
    originals_by_path = {item["relative_path"]: item for item in inputs}
    require({row["source_path"] for row in inventory} == set(originals_by_path),
            "Inventory source paths differ from loaded inputs")
    for item in inventory:
        original = originals_by_path[item["source_path"]]
        require(item["source_sha256"] == original["sha256"], "Inventory source checksum mismatch")
        if "repository_path" in item:
            require(item["repository_path"] == original["repository_path"],
                    "Inventory physical repository path is stale after relocation")
    require(len({cell["cell_id"] for cell in cells}) == len(cells), "Duplicate exported cell IDs")
    observed = set()
    counts = Counter()
    value_by_id = {value["mapping_id"]: value for value in values}
    require(len(value_by_id) == len(values), "Duplicate value mapping IDs")
    for cell in cells:
        location = (cell["source_path"], int(cell["source_record"]), int(cell["column_index"]))
        require(location in source and location not in observed, f"Invalid/duplicate source coordinate: {location}")
        observed.add(location)
        original = source[location]
        require(cell["raw_value"] == original["value"], f"Raw cell changed: {location}")
        require(cell["source_column"] == original["header"], f"Original header changed: {location}")
        require(cell["source_sha256"] == original["sha256"], f"Source checksum mismatch: {location}")
        require(int(cell["source_line_start"]) == original["line_start"] and
                int(cell["source_line_end"]) == original["line_end"], f"Physical CSV line mismatch: {location}")
        require(cell["mapping_id"] in value_by_id, f"Missing value mapping: {location}")
        mapping = value_by_id[cell["mapping_id"]]
        for key in ("raw_value", "status", "value_kind", "term_id", "unit_id", "numeric_value", "normalized_value", "assertions"):
            require(cell.get(key, "") == mapping.get(key, ""), f"Cell/value disagreement for {key}: {location}")
        check_terms(cell, str(location))
        counts[cell["mapping_id"]] += 1
    require(observed == set(source), "Not every original source coordinate is exported")
    for value in values:
        require(counts[value["mapping_id"]] == int(value["cell_count"]), "Mapping expansion count mismatch")
        check_terms(value, value["mapping_id"])
    for column in columns:
        check_terms(column, column["source_path"] + "/" + column["source_column"], column=True)
    require(dict(Counter(cell["status"] for cell in cells)) == report["cell_status_counts"], "QA status totals differ from saved cells")
    code_dir = Path(__file__).parent
    for name, expected in report["rule_file_sha256"].items():
        require(digest((code_dir / name).read_bytes()) == expected, f"Rules changed since export: {name}")
    for name, expected in report["catalog_file_sha256"].items():
        require(digest((catalog.DATA / name).read_bytes()) == expected, f"Vocabulary changed since export: {name}")
    for name, expected in report.get("evidence_file_sha256", {}).items():
        require(digest((code_dir.parent / name).read_bytes()) == expected, f"Evidence changed since export: {name}")
    project_exports_checked = validate_project_exports(
        output, report, {"cell_annotations": cells, "value_mappings": values,
                         "column_mappings": columns, "row_annotations": rows,
                         "inventory": inventory, "curation_queue": curation}, project_output_root)
    result = {"validation": "passed", "source_cells_checked": len(cells),
              "source_records_checked": len(rows), "column_instances_checked": len(columns),
              "value_mappings_checked": len(values),
              "project_exports_checked": project_exports_checked,
              "validator_sha256": digest(Path(__file__).read_bytes()),
              "checks": ["TSV round-trip raw values and multiline headers", "complete unique cell coordinates",
                         "original source hashes and physical lines", "saved value dictionary expansion",
                         "all assigned ontology IDs current and verified", "review candidates not assigned",
                         "rule and vocabulary hashes match the export",
                         "project exports equal their aggregate subsets"]}
    (output / "export_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--project-output-root", type=Path,
                        help="override per-project artifact root; defaults to paths in aggregate QA metadata")
    args = parser.parse_args()
    print(json.dumps(validate(args.output, args.data_root, args.project_output_root), indent=2))


if __name__ == "__main__":
    main()
