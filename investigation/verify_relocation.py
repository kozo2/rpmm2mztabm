#!/usr/bin/env python3
"""Check moved bytes, stable scientific annotations and project-local exports.

The checked-in manifest contains pre-move hashes and semantic invariants. No
temporary baseline directory, network access or optional package is required.
"""

import csv
import argparse
import hashlib
import json
from pathlib import Path

from phenodata_ontology.validate_export import validate

REPO = Path(__file__).resolve().parent.parent
REPORT = REPO / "reports" / "relocation"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuilt-tables-root", type=Path,
                        help="check a fresh export_tables --output-root result against every original table hash")
    args = parser.parse_args()
    manifest = json.loads((REPORT / "manifest.json").read_text())
    counts = {}
    for entry in manifest["files"]:
        path = REPO / entry["new_path"]
        require(path.is_file(), "Missing moved file: " + str(path))
        mode = entry["validation_mode"]
        counts[mode] = counts.get(mode, 0) + 1
        if mode == "byte_identity":
            require(hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"],
                    "Moved bytes changed: " + entry["new_path"])
        elif mode == "table_report_paths":
            actual = rows(path)
            for record in actual:
                for key in ("toml", "parquet"):
                    target = REPO / record.pop(key)
                    require(target.is_file(), "Missing table reference: " + str(target))
                    require(target.parent == REPO / "projects" / record["project"] / "tables",
                            "Table reference points outside its project")
            require(actual == entry["expected_rows"], "Table report scientific content changed")
        elif mode == "inventory_with_physical_paths":
            actual = rows(path)
            for record in actual:
                target = REPO / record.pop("repository_path")
                project, tail = record["source_path"].split("/", 1)
                require(target == REPO / "projects" / project / "data" / tail,
                        "Logical-to-physical source resolution changed")
                require(hashlib.sha256(target.read_bytes()).hexdigest() == record["source_sha256"],
                        "Annotation source content changed")
            require(actual == entry["expected_rows"], "Inventory scientific content changed")
        elif mode in {"qa_invariants", "export_validation_metrics"}:
            actual = json.loads(path.read_text())
            for key, value in entry["expected_metrics"].items():
                require(actual.get(key) == value, "Reproduction metric changed: " + key)
        else:
            raise ValueError("Unknown relocation validation mode: " + mode)
    annotation = validate()
    if args.rebuilt_tables_root:
        tables = [entry for entry in manifest["files"]
                  if entry["new_path"].startswith("projects/") and "/tables/" in entry["new_path"]]
        for entry in tables:
            rebuilt = args.rebuilt_tables_root / Path(entry["new_path"]).relative_to("projects")
            require(rebuilt.is_file(), "Rebuilt table missing: " + str(rebuilt))
            require(hashlib.sha256(rebuilt.read_bytes()).hexdigest() == entry["sha256"],
                    "Regenerated table bytes differ: " + str(rebuilt))
        table_reproduction = {"compared_files": len(tables), "byte_identical": len(tables),
                              "semantic_equal": len(tables), "differences": [],
                              "generated_on": "2026-10-27",
                              "generator_sha256": hashlib.sha256((REPO / "investigation/export_tables.py").read_bytes()).hexdigest(),
                              "requirements_sha256": hashlib.sha256((REPO / "investigation/requirements-export.txt").read_bytes()).hexdigest()}
        (REPORT / "table_reproduction.json").write_text(json.dumps(table_reproduction, indent=2) + "\n")
    table_reproduction = json.loads((REPORT / "table_reproduction.json").read_text())
    require(table_reproduction["byte_identical"] == table_reproduction["compared_files"] == 212,
            "Table regeneration was not byte-identical")
    require(table_reproduction["generator_sha256"] == hashlib.sha256(
        (REPO / "investigation/export_tables.py").read_bytes()).hexdigest(), "Table generator changed since reproduction check")
    require(table_reproduction["requirements_sha256"] == hashlib.sha256(
        (REPO / "investigation/requirements-export.txt").read_bytes()).hexdigest(), "Table dependencies changed since reproduction check")
    result = {"validation": "passed", "moved_files": len(manifest["files"]),
              "preservation_checks": counts, "project_folders": len(list((REPO / "projects").glob("RPMM*"))),
              "annotation_cells_checked": annotation["source_cells_checked"],
              "project_annotation_exports_checked": annotation["project_exports_checked"],
              "regenerated_table_files_byte_identical": table_reproduction["byte_identical"],
              "logical_cell_and_mapping_ids": "unchanged (complete annotation TSV byte identity)",
              "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (REPORT / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
