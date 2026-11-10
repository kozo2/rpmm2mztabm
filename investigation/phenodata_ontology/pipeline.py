"""Deterministic, auditable column/value/cell annotation export (stdlib only)."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import os
import sys

from . import catalog
from .columns import classify_column
from .values import annotate_value
from .provenance import load_discovery, manifest_records, refresh_discovery

INVESTIGATION = Path(__file__).resolve().parents[1]
REPOSITORY = INVESTIGATION.parent
DEFAULT_DATA = REPOSITORY / "projects"
DEFAULT_OUTPUT = REPOSITORY / "reports" / "phenodata"
DEFAULT_PROJECT_OUTPUT = DEFAULT_DATA


def packed(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def digest(value):
    raw = value if isinstance(value, bytes) else packed(value).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def load_inputs(data_root, projects=None):
    """Read project/data trees while retaining stable RPMM/... logical keys.

    Flat RPMM/... data roots remain readable for archived inputs and fixtures.
    Generated project artifacts are never searched when a data subfolder exists.
    """
    data_root = Path(data_root)
    inputs = []
    paths = []
    for folder in sorted(data_root.glob("RPMM*")):
        if not folder.is_dir() or (projects and folder.name not in projects):
            continue
        source_root = folder / "data" if (folder / "data").is_dir() else folder
        paths.extend((path, folder.name, source_root) for path in source_root.rglob("*")
                     if path.is_file() and path.name.lower() == "phenodata.csv")
    for path, project, source_root in sorted(paths, key=lambda entry: str(entry[0])):
        if not path.is_file() or path.name.lower() != "phenodata.csv":
            continue
        relative = project + "/" + path.relative_to(source_root).as_posix()
        if not project.startswith("RPMM") or (projects and project not in projects):
            continue
        raw = path.read_bytes()
        records = []
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            headers = next(reader)
            if len(headers) != len(set(headers)):
                raise ValueError(f"Duplicate CSV headers are ambiguous: {relative}")
            previous_line = reader.line_num
            for record_number, cells in enumerate(reader, 2):
                if len(cells) != len(headers):
                    raise ValueError(f"CSV width mismatch: {relative} record {record_number}")
                records.append({"source_record": record_number,
                                "source_line_start": previous_line + 1,
                                "source_line_end": reader.line_num,
                                "cells": cells})
                previous_line = reader.line_num
        inputs.append({"path": path, "relative_path": relative, "project": project,
                       "source_root": source_root,
                       "repository_path": (path.relative_to(REPOSITORY).as_posix()
                                           if REPOSITORY in path.parents else str(path.resolve())),
                       "sha256": digest(raw), "md5": hashlib.md5(raw).hexdigest(),
                       "headers": headers, "records": records})
    if not inputs:
        raise ValueError(f"No phenodata.csv files found under {data_root}")
    return inputs


def classify_row(row):
    """Identify explicit standards/references; do not infer species from project."""
    values = [value.strip() for value in row.values()]
    if not any(values):
        return "empty", "all source cells are empty"
    lowers = {value.lower() for value in values}
    if ({"standards mixture", "standard mixture", "stds"} & lowers
            or row.get("sample", "").strip().lower() == "standard"
            or row.get("type", "").strip().lower() in {"sd", "alkanemix"}):
        return "standard_material", "explicit STDs/Standards mixture/standard type marker"
    if any("freeze dry col-0" in value.lower() for value in values):
        return "reference_material", "explicit Freeze Dry Col-0 reference marker"
    if "fd" in lowers and sum(value.lower() == "fd" for value in values) >= 3:
        return "reference_material", "FD marker in at least three metadata fields"
    return "experimental_record", "no explicit standard/reference marker; biological identity not inferred"


def pick_identifier(headers, cells):
    priorities = ["gcid", "sample id", "sampleid", "numbering", "filename", "lcid",
                  "msp", "fddap", "sample name", "samplename"]
    for name in priorities:
        for header, value in zip(headers, cells):
            if header.strip().lower() == name and value.strip():
                return header, value
    return "", ""


def validate_annotation(annotation, force_review=False):
    """Only emit current, source-verified ontology identifiers as assigned IDs."""
    result = dict(annotation)
    result.setdefault("assertions", [])
    result.setdefault("status", "accepted" if result.get("term_id") or result.get("unit_id") else "local")
    if force_review:
        result["status"] = "needs_review"
    problems = []
    for key in ("term_id", "unit_id"):
        term_id = result.get(key, "")
        if not term_id:
            continue
        term = catalog.get(term_id)
        if not term or not catalog.is_current(term_id):
            result["candidate_" + key] = term_id
            result[key] = ""
            problems.append(f"{key} {term_id} is not in the verified current catalog")
        else:
            stem = "term" if key == "term_id" else "unit"
            result[stem + "_label"] = term.get("label", "")
            result[stem + "_iri"] = term.get("iri", "")
            result[stem + "_source_id"] = term.get("source_id", "")
            result[stem + "_ontology_version"] = term.get("ontology_version", "")
    checked = []
    for assertion in result["assertions"]:
        if not isinstance(assertion, dict):
            problems.append("invalid assertion object")
            continue
        validated = validate_annotation(assertion,
                                        force_review=result.get("status") == "needs_review")
        if validated.get("status") == "needs_review":
            problems.append("a component assertion needs review")
        checked.append(validated)
    result["assertions"] = checked
    if problems:
        result["status"] = "needs_review"
        result["reason"] = "; ".join(filter(None, [result.get("reason", ""), *problems]))
    if result.get("status") == "needs_review":
        # A later invalid component can downgrade the parent after an earlier
        # component was verified. Do not leak that earlier ID as an assignment.
        result["assertions"] = [validate_annotation(assertion, force_review=True)
                                for assertion in result["assertions"]]
        for key in ("term_id", "unit_id"):
            if result.get(key):
                result["candidate_" + key] = result[key]
                result[key] = ""
                for suffix in ("label", "iri", "source_id", "ontology_version"):
                    stem = "term" if key == "term_id" else "unit"
                    result.pop(stem + "_" + suffix, None)
    return result


def validate_column(spec):
    result = dict(spec)
    for key in ("column_term_id", "unit_id"):
        term_id = result.get(key, "")
        term = catalog.get(term_id) if term_id else None
        if term_id and (not term or not catalog.is_current(term_id)):
            result["candidate_" + key] = term_id
            result[key] = ""
            result["status"] = "needs_review"
            result["reason"] = (result.get("reason", "") +
                                f"; {term_id} requires current catalog verification").strip("; ")
        elif term:
            prefix = "column_term" if key == "column_term_id" else "unit"
            result[prefix + "_label"] = term.get("label", "")
            result[prefix + "_iri"] = term.get("iri", "")
            result[prefix + "_source_id"] = term.get("source_id", "")
            result[prefix + "_ontology_version"] = term.get("ontology_version", "")
    return result


@lru_cache(maxsize=1)
def _vocabulary():
    return catalog.OntologyCatalog()


@lru_cache(maxsize=None)
def _candidate_lookup(spelling, prefix):
    return _vocabulary().find_exact(spelling, prefix)


def add_candidates(annotation, spec, raw_value):
    """Use current native labels/EXACT synonyms as suggestions, never approvals."""
    result = dict(annotation)
    if result.get("status") != "needs_review":
        return result
    prefix = {"organism_part": "PO", "organism": "NCBITaxon",
              "treatment": "PECO", "phenotype": "PATO"}.get(spec.get("semantic_role"))
    if not prefix:
        return result
    terms = {}
    for spelling in {raw_value.strip(), result.get("normalized_value", "")}:
        for term in _candidate_lookup(spelling, prefix):
            terms[term["id"]] = {"term_id": term["id"], "label": term["label"],
                                 "iri": term["iri"], "source_id": term["source_id"],
                                 "match": "native label or EXACT synonym", "status": "candidate"}
    if terms:
        result["candidate_terms"] = sorted(terms.values(), key=lambda term: term["term_id"])
    return result


def write_tsv(path, rows, preferred=()):
    keys = list(preferred)
    for key in sorted({key for row in rows for key in row}):
        if key not in keys:
            keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: packed(value) if isinstance(value, (dict, list)) else
                             ("" if value is None else value) for key, value in row.items()})
    return keys


def base_evidence(item):
    return {"project": item["project"], "source_path": item["relative_path"],
            "source_sha256": item["sha256"]}


def detect_conflicts(item, source_record, row):
    issues = []
    if item["project"] == "RPMM0017" and "valid02_gc" in item["relative_path"]:
        import re
        note = row.get("note", "")
        match = re.search(r"leaf\s*:\s*fruit\s*=\s*(\d+)\s*:\s*(\d+)", note, re.I)
        if match:
            try:
                mismatch = (float(match[1]) != float(row["mixtureratioleaf"]) or
                            float(match[2]) != float(row["mixtureratiofruit"]))
            except (ValueError, KeyError):
                mismatch = True
            if mismatch:
                issues.append({"issue_code": "conflicting_mixture_ratio", "severity": "error",
                               "columns": ["mixtureratiofruit", "mixtureratioleaf", "note"],
                               "reason": "leaf:fruit ratio in note disagrees with named numeric columns",
                               "suggested_action": "Consult original protocol/publication; do not choose either ratio automatically"})
    return issues


def run(data_root=DEFAULT_DATA, output=DEFAULT_OUTPUT, projects=None, project_output_root=None):
    data_root, output = Path(data_root).resolve(), Path(output).resolve()
    inputs = load_inputs(data_root, projects)
    project_output_root = Path(project_output_root).resolve() if project_output_root else (
        data_root if output == DEFAULT_OUTPUT.resolve() else output / "projects")
    destinations = [output] + [project_output_root / item["project"] / "phenodata" for item in inputs]
    for destination in destinations:
        if destination == data_root or any(item["source_root"] == destination or
                item["source_root"] in destination.parents for item in inputs):
            raise ValueError("Annotation output must be outside the immutable project data folders")
    output.mkdir(parents=True, exist_ok=True)
    columns, cells, row_records, inventory = [], [], [], []
    mappings = {}
    curation = []
    manifest = manifest_records(INVESTIGATION)
    discovery = load_discovery()
    if discovery.get("status") != "sample_links_verified":
        curation.append({"scope": "external_source", "severity": "warning",
                         "issue_code": "sample_level_rdf_context_unavailable",
                         "reason": discovery.get("limitation", "No sample-level RDF annotations have been verified"),
                         "suggested_action": "Retry PMM discovery; verify experiment/file/sample joins before reusing RDF annotations"})
    for item in inputs:
        base = base_evidence(item)
        specs = [validate_column(classify_column(item["project"], item["relative_path"], header))
                 for header in item["headers"]]
        for index, (header, spec) in enumerate(zip(item["headers"], specs), 1):
            record = {**base, "column_index": index, "source_column": header, **spec}
            columns.append(record)
            if spec.get("status") == "needs_review":
                curation.append({**base, "scope": "column", "column_index": index,
                                 "source_column": header, "severity": "warning",
                                 "issue_code": "column_semantics", "reason": spec.get("reason", ""),
                                 "suggested_action": "Confirm column definition/unit against source protocol"})
        row_kinds = Counter()
        for record in item["records"]:
            row = dict(zip(item["headers"], record["cells"]))
            row_kind, row_evidence = classify_row(row)
            row_kinds[row_kind] += 1
            id_column, sample_id = pick_identifier(item["headers"], record["cells"])
            location = {key: value for key, value in record.items() if key != "cells"}
            row_id = digest([base["source_path"], base["source_sha256"], record["source_record"]])
            row_records.append({**base, **location, "row_id": row_id, "row_kind": row_kind,
                                "row_kind_evidence": row_evidence, "sample_id_column": id_column,
                                "sample_id": sample_id})
            conflicts = detect_conflicts(item, record["source_record"], row)
            for conflict in conflicts:
                curation.append({**base, **location, "scope": "row", "sample_id": sample_id,
                                 **conflict})
            conflict_columns = {name for conflict in conflicts for name in conflict["columns"]}
            for index, (header, raw_value, spec) in enumerate(zip(item["headers"], record["cells"], specs), 1):
                annotated = annotate_value(item["project"], item["relative_path"], header,
                                           raw_value, row, spec)
                if (spec.get("status") == "needs_review" and
                        annotated.get("status") not in {"missing", "not_applicable"}):
                    annotated["status"] = "needs_review"
                    annotated["reason"] = "; ".join(filter(None, [annotated.get("reason", ""),
                        "Column definition requires review: " + spec.get("reason", "")]))
                if row_kind == "empty":
                    annotated = {"value_kind": "missing", "status": "missing",
                                 "reason": "empty source row retained", "rule_id": "empty_row"}
                elif header in conflict_columns:
                    annotated = {**annotated, "status": "needs_review",
                                 "reason": "conflicting_mixture_ratio: numeric columns and note disagree"}
                annotated = validate_annotation(annotated)
                annotated = add_candidates(annotated, spec, raw_value)
                annotated.setdefault("value_kind", spec.get("value_kind", "unresolved"))
                annotated.setdefault("status", "needs_review")
                annotated.setdefault("normalized_value", raw_value.strip())
                # Output-dependent key prevents context-sensitive values from sharing incompatible mappings.
                mapping_id = digest([base["source_path"], index, raw_value, row_kind, annotated])
                if mapping_id not in mappings:
                    mappings[mapping_id] = {**base, "mapping_id": mapping_id,
                                           "column_index": index, "source_column": header,
                                           "standard_name": spec.get("standard_name", ""),
                                           "semantic_role": spec.get("semantic_role", ""),
                                           "raw_value": raw_value, "row_kind": row_kind,
                                           **annotated, "cell_count": 0, "source_records": []}
                mappings[mapping_id]["cell_count"] += 1
                mappings[mapping_id]["source_records"].append(record["source_record"])
                cells.append({**base, **location, "row_id": row_id,
                              "cell_id": digest([row_id, index]), "column_index": index,
                              "source_column": header, "raw_value": raw_value,
                              "standard_name": spec.get("standard_name", ""),
                              "semantic_role": spec.get("semantic_role", ""),
                              "row_kind": row_kind, "sample_id": sample_id,
                              "mapping_id": mapping_id, **annotated})
        source_manifest = manifest.get(item["relative_path"], {})
        if source_manifest and source_manifest.get("md5") != item["md5"]:
            curation.append({**base, "scope": "source_file", "severity": "warning",
                             "issue_code": "source_manifest_checksum_mismatch",
                             "reason": "Current file MD5 differs from published manifest MD5; annotations identify the actual source by SHA256",
                             "published_md5": source_manifest.get("md5", ""), "actual_md5": item["md5"],
                             "download_time_md5_match": source_manifest.get("md5_match", ""),
                             "suggested_action": "Resolve the pre-existing download checksum discrepancy before asserting identity with remote sample metadata"})
        inventory.append({**base, "columns": len(item["headers"]),
                          "repository_path": item["repository_path"],
                          "data_records": len(item["records"]), "row_kind_counts": dict(row_kinds),
                          "download_url": source_manifest.get("downloadURL", ""),
                          "source_file_iri": source_manifest.get("file_iri", ""),
                          "manifest_md5": source_manifest.get("md5", ""), "current_md5": item["md5"],
                          "manifest_md5_matches": item["md5"] == source_manifest.get("md5", "")})
    values = list(mappings.values())
    for value in values:
        if value["status"] == "needs_review":
            curation.append({key: value[key] for key in (
                "project", "source_path", "source_sha256", "mapping_id", "column_index",
                "source_column", "raw_value", "cell_count", "source_records", "reason")})
            curation[-1].update({"scope": "value", "severity": "warning", "issue_code": "value_semantics",
                                 "suggested_action": "Review source context and candidate term definitions"})
            for key in ("candidate_term_id", "candidate_unit_id", "candidate_terms", "assertions"):
                if value.get(key):
                    curation[-1][key] = value[key]
    qa = verify(inputs, columns, values, cells, row_records)
    qa.update({"generated_at_utc": datetime.now(timezone.utc).isoformat(),
               "python_version": sys.version.split()[0],
               "scope_projects": sorted(projects) if projects else "all matching RPMM projects",
               "rule_file_sha256": {path.name: digest(path.read_bytes()) for path in
                   [Path(__file__), Path(__file__).with_name("columns.py"), Path(__file__).with_name("values.py"),
                    Path(__file__).with_name("catalog.py"), Path(__file__).with_name("provenance.py")]},
               "catalog_file_sha256": {path.name: digest(path.read_bytes()) for path in
                   [catalog.DATA / "ontology_terms.json", catalog.DATA / "ontology_sources.json"]},
               "evidence_file_sha256": {str(path.relative_to(INVESTIGATION)): digest(path.read_bytes())
                   for path in [INVESTIGATION / "projects.json", INVESTIGATION / "manifest.csv",
                                catalog.DATA / "rdf_discovery.json"] if path.exists()},
               "rdf_discovery_status": discovery.get("status"),
               "rdf_discovery_limitation": discovery.get("limitation", ""),
               "cell_status_counts": dict(Counter(row["status"] for row in cells)),
               "cell_kind_counts": dict(Counter(row["value_kind"] for row in cells)),
               "column_status_counts": dict(Counter(row["status"] for row in columns)),
               "value_status_counts": dict(Counter(row["status"] for row in values)),
               "row_kind_counts": dict(Counter(row["row_kind"] for row in row_records)),
               "curation_items": len(curation),
               "cells_with_assigned_term": sum(bool(row.get("term_id")) for row in cells),
               "cells_with_assigned_unit": sum(bool(row.get("unit_id")) for row in cells),
               "columns_with_assigned_term": sum(bool(row.get("column_term_id")) for row in columns)})
    tables = dict([("inventory", inventory), ("column_mappings", columns),
                       ("value_mappings", values), ("cell_annotations", cells),
                       ("row_annotations", row_records), ("curation_queue", curation)])
    fields = {}
    for name, rows in tables.items():
        fields[name] = write_tsv(output / (name + ".tsv"), rows,
                                ("project", "source_path", "source_sha256"))
    coverage = []
    for project in sorted({item["project"] for item in inputs}):
        project_cells = [cell for cell in cells if cell["project"] == project]
        project_values = [value for value in values if value["project"] == project]
        counts = Counter(cell["status"] for cell in project_cells)
        coverage.append({"project": project,
                         "files": sum(item["project"] == project for item in inputs),
                         "cells": len(project_cells), "unique_value_mappings": len(project_values),
                         "cells_with_term_id": sum(bool(cell.get("term_id")) for cell in project_cells),
                         "cells_with_unit_id": sum(bool(cell.get("unit_id")) for cell in project_cells),
                         **{status: counts[status] for status in ["accepted", "local", "needs_review", "missing", "not_applicable"]}})
    write_tsv(output / "coverage_by_project.tsv", coverage, ("project", "files", "cells"))
    role_coverage = []
    for role in sorted({cell["semantic_role"] for cell in cells}):
        role_cells = [cell for cell in cells if cell["semantic_role"] == role]
        role_values = [value for value in values if value["semantic_role"] == role]
        role_coverage.append({"semantic_role": role, "cells": len(role_cells), "unique_value_mappings": len(role_values),
                              "cells_with_assigned_term": sum(bool(cell.get("term_id")) for cell in role_cells),
                              "values_with_assigned_term": sum(bool(value.get("term_id")) for value in role_values),
                              "cells_with_assigned_unit": sum(bool(cell.get("unit_id")) for cell in role_cells),
                              "values_with_assigned_unit": sum(bool(value.get("unit_id")) for value in role_values),
                              "cells_needing_review": sum(cell["status"] == "needs_review" for cell in role_cells)})
    write_tsv(output / "coverage_by_role.tsv", role_coverage, ("semantic_role", "cells", "unique_value_mappings"))
    qa["project_exports"] = {}
    for project in sorted({item["project"] for item in inputs}):
        folder = project_output_root / project / "phenodata"
        folder.mkdir(parents=True, exist_ok=True)
        qa["project_exports"][project] = {"directory": os.path.relpath(folder, output)}
        for name, table in tables.items():
            write_tsv(folder / (name + ".tsv"), [row for row in table if row.get("project") == project], fields[name])
        write_tsv(folder / "coverage_by_project.tsv", [row for row in coverage if row["project"] == project],
                  ("project", "files", "cells"))
        selected_cells = [cell for cell in cells if cell["project"] == project]
        selected_values = [value for value in values if value["project"] == project]
        per_role = []
        for role in sorted({cell["semantic_role"] for cell in selected_cells}):
            cs = [cell for cell in selected_cells if cell["semantic_role"] == role]
            vs = [value for value in selected_values if value["semantic_role"] == role]
            per_role.append({"semantic_role": role, "cells": len(cs), "unique_value_mappings": len(vs),
                             "cells_with_assigned_term": sum(bool(cell.get("term_id")) for cell in cs),
                             "values_with_assigned_term": sum(bool(value.get("term_id")) for value in vs),
                             "cells_with_assigned_unit": sum(bool(cell.get("unit_id")) for cell in cs),
                             "values_with_assigned_unit": sum(bool(value.get("unit_id")) for value in vs),
                             "cells_needing_review": sum(cell["status"] == "needs_review" for cell in cs)})
        write_tsv(folder / "coverage_by_role.tsv", per_role,
                  ("semantic_role", "cells", "unique_value_mappings"))
        provenance = {"project": project, "aggregate_report": os.path.relpath(output / "qa_report.json", folder),
                      "files": sum(item["project"] == project for item in inputs),
                      "source_cells": len(selected_cells),
                      "source_records": sum(row["project"] == project for row in row_records),
                      "source_column_instances": sum(column["project"] == project for column in columns),
                      "unique_value_mappings": len(selected_values),
                      **{key: qa[key] for key in ("rule_file_sha256", "catalog_file_sha256", "evidence_file_sha256")}}
        (folder / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n")
    (output / "qa_report.json").write_text(json.dumps(qa, indent=2, ensure_ascii=False) + "\n")
    return qa


def verify(inputs, columns, values, cells, row_records):
    expected_cells = sum(len(item["headers"]) * len(item["records"]) for item in inputs)
    expected_rows = sum(len(item["records"]) for item in inputs)
    expected_columns = sum(len(item["headers"]) for item in inputs)
    assert len(cells) == expected_cells, "source cell coverage mismatch"
    assert len(row_records) == expected_rows, "source row coverage mismatch"
    assert len(columns) == expected_columns, "source column coverage mismatch"
    assert len({cell["cell_id"] for cell in cells}) == expected_cells, "duplicate cell identifiers"
    assert sum(value["cell_count"] for value in values) == expected_cells, "value expansion mismatch"
    assert {cell["mapping_id"] for cell in cells} == {value["mapping_id"] for value in values}
    originals = {(item["relative_path"], record["source_record"], index): raw
                 for item in inputs for record in item["records"]
                 for index, raw in enumerate(record["cells"], 1)}
    for cell in cells:
        assert cell["raw_value"] == originals[(cell["source_path"], cell["source_record"], cell["column_index"])], "raw value changed"
        assert cell["status"] in {"accepted", "local", "needs_review", "missing", "not_applicable"}
    for item in inputs:
        assert digest(item["path"].read_bytes()) == item["sha256"], "input file modified during run"
    return {"validation": "passed", "files": len(inputs),
            "projects": len({item["project"] for item in inputs}),
            "source_records": expected_rows,
            "nonempty_records": sum(row["row_kind"] != "empty" for row in row_records),
            "source_column_instances": expected_columns,
            "distinct_raw_headers": len({header for item in inputs for header in item["headers"]}),
            "source_cells": expected_cells, "unique_value_mappings": len(values),
            "checks": ["every source record/column/cell represented", "unique cell identity",
                       "raw values preserved", "mapping expansion is complete", "input hashes unchanged"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("projects", nargs="*", help="optional RPMM project IDs")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--project-output-root", type=Path,
                        help="root for RPMM####/phenodata exports; custom aggregate outputs default to OUTPUT/projects")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--refresh-rdf", action="store_true", help="retry public PMM file-link discovery (network required)")
    args = parser.parse_args()
    if args.projects and args.output == DEFAULT_OUTPUT:
        parser.error("Project subsets require --output so a pilot cannot replace the full export")
    if args.refresh_rdf:
        inputs = load_inputs(args.data_root.resolve(), set(args.projects) or None)
        manifest = manifest_records(INVESTIGATION)
        refresh_discovery([manifest[item["relative_path"]]["file_iri"] for item in inputs
                           if item["relative_path"] in manifest])
    report = run(args.data_root.resolve(), args.output.resolve(), set(args.projects) or None,
                 args.project_output_root)
    print(json.dumps(report, indent=2, ensure_ascii=False))
