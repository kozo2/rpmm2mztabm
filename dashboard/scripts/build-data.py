#!/usr/bin/env python3
"""Build a deterministic, provenance-preserving RPMM dashboard data snapshot.

The column mapping TSV defines the set of available source columns. Values are
joined by project, logical source path and one-based column index, never by a
header or ontology label alone. No annotations are inferred by this generator.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import shutil


DASHBOARD = Path(__file__).resolve().parents[1]
REPOSITORY = DASHBOARD.parent
DEFAULT_REPORTS = REPOSITORY / "reports" / "phenodata"
DEFAULT_TITLES = REPOSITORY / "investigation" / "projects.json"
DEFAULT_OUTPUT = DASHBOARD / "snapshot"


def stable_json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def integer(value, field):
    try:
        result = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid integer in {field}: {value!r}") from error
    if result < 0:
        raise ValueError(f"Negative count in {field}: {result}")
    return result


def column_key(row):
    index = integer(row["column_index"], "column_index")
    if index == 0:
        raise ValueError("Source column indexes must be one-based")
    return "|".join((row["project"], row["source_path"], str(index)))


def term(row, prefix):
    identifier = row.get(prefix + "_id", "")
    if not identifier:
        return None
    return {
        "id": identifier,
        "label": row.get(prefix + "_label", ""),
        "iri": row.get(prefix + "_iri", ""),
        "version": row.get(prefix + "_ontology_version", ""),
        "sourceId": row.get(prefix + "_source_id", ""),
    }


def json_array(value, name):
    parsed = json.loads(value) if value else []
    if not isinstance(parsed, list):
        raise ValueError(f"{name} must contain a JSON array")
    return parsed


def title_text(value):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]*>", "", value))).strip()


def source_label(path):
    """Use portable provenance paths even when a caller supplies absolute roots."""
    path = Path(path).resolve()
    try:
        return path.relative_to(REPOSITORY).as_posix()
    except ValueError:
        return path.name


def build_snapshot(reports_root=DEFAULT_REPORTS, titles_path=DEFAULT_TITLES):
    reports_root, titles_path = Path(reports_root), Path(titles_path)
    column_path = reports_root / "column_mappings.tsv"
    value_path = reports_root / "value_mappings.tsv"
    inventory_path = reports_root / "inventory.tsv"
    columns_raw = read_tsv(column_path)
    values_raw = read_tsv(value_path)
    inventory_raw = read_tsv(inventory_path)
    titles_raw = json.loads(titles_path.read_text(encoding="utf-8"))
    titles = {record["p"].rsplit("/", 1)[-1]: title_text(record.get("label", ""))
              for record in titles_raw}

    source_files = [{"path": source_label(path), "sha256": sha256(path)}
                    for path in (column_path, value_path, inventory_path, titles_path)]
    source_files.sort(key=lambda entry: entry["path"])
    dataset_sha = hashlib.sha256(stable_json(source_files).encode("utf-8")).hexdigest()

    source_by_path = {}
    files = []
    for raw in inventory_raw:
        key = (raw["project"], raw["source_path"])
        if key in source_by_path:
            raise ValueError(f"Duplicate inventory source: {key}")
        record = {
            "project": raw["project"], "sourcePath": raw["source_path"],
            "sourceSha256": raw["source_sha256"],
            "repositoryPath": raw.get("repository_path", ""),
            "dataRecords": integer(raw["data_records"], "data_records"),
            "columnCount": integer(raw["columns"], "inventory columns"),
            "downloadUrl": raw.get("download_url", ""),
            "sourceFileIri": raw.get("source_file_iri", ""),
            "manifestMd5Matches": raw.get("manifest_md5_matches", "").lower() == "true",
        }
        source_by_path[key] = record
        files.append(record)

    columns = []
    column_by_key = {}
    for raw in columns_raw:
        key = column_key(raw)
        if key in column_by_key:
            raise ValueError(f"Duplicate source column: {key}")
        source = source_by_path.get((raw["project"], raw["source_path"]))
        if source is None or source["sourceSha256"] != raw["source_sha256"]:
            raise ValueError(f"Column has no matching inventory source/hash: {key}")
        concept = term(raw, "column_term")
        record = {
            "key": key, "project": raw["project"], "sourcePath": raw["source_path"],
            "sourceSha256": raw["source_sha256"], "columnIndex": integer(raw["column_index"], "column_index"),
            "sourceColumn": raw["source_column"], "standardName": raw["standard_name"],
            "semanticRole": raw["semantic_role"], "valueKind": raw["value_kind"],
            "status": raw["status"], "columnTerm": concept, "unit": term(raw, "unit"),
            "groupKey": "term:" + concept["id"] if concept else "field:" + raw["standard_name"],
            "referenceEvent": raw.get("reference_event", ""),
            "reason": raw.get("reason", ""), "evidence": raw.get("evidence", ""),
            "ruleId": raw.get("rule_id", ""),
        }
        columns.append(record)
        column_by_key[key] = record

    values_by_column = defaultdict(list)
    seen_mapping_ids = set()
    for raw in values_raw:
        key = column_key(raw)
        column = column_by_key.get(key)
        if column is None:
            raise ValueError(f"Value refers to a nonexistent source column: {key}")
        if (raw["source_column"] != column["sourceColumn"] or
                raw["source_sha256"] != column["sourceSha256"]):
            raise ValueError(f"Value and column provenance disagree: {key}")
        mapping_id = raw["mapping_id"]
        if not mapping_id or mapping_id in seen_mapping_ids:
            raise ValueError(f"Missing or duplicate value mapping ID: {mapping_id}")
        seen_mapping_ids.add(mapping_id)
        count = integer(raw["cell_count"], "cell_count")
        records = json_array(raw.get("source_records", ""), "source_records")
        if len(records) != count or any(not isinstance(record, int) or record < 2 for record in records):
            raise ValueError(f"Source records do not match cell_count: {mapping_id}")
        row = {
            "mappingId": mapping_id, "rawValue": raw["raw_value"],
            "normalizedValue": raw.get("normalized_value", ""),
            # Decimal strings preserve the precision of the recorded annotation.
            "numericValue": raw.get("numeric_value", ""),
            "status": raw["status"], "valueKind": raw["value_kind"],
            "term": term(raw, "term"), "unit": term(raw, "unit"),
            "candidateTermId": raw.get("candidate_term_id", ""),
            "candidateUnitId": raw.get("candidate_unit_id", ""),
            "candidateTerms": json_array(raw.get("candidate_terms", ""), "candidate_terms"),
            "assertions": json_array(raw.get("assertions", ""), "assertions"),
            "referenceEvent": raw.get("reference_event", ""), "cellCount": count,
            "rowKind": raw.get("row_kind", ""), "sourceRecords": records,
            "reason": raw.get("reason", ""), "evidence": raw.get("evidence", ""),
            "ruleId": raw.get("rule_id", ""),
        }
        if row["status"] == "needs_review" and (row["term"] or row["unit"]):
            raise ValueError(f"A review candidate is incorrectly assigned: {mapping_id}")
        values_by_column[key].append(row)

    for column in columns:
        values = values_by_column[column["key"]]
        source = source_by_path[(column["project"], column["sourcePath"])]
        if sum(value["cellCount"] for value in values) != source["dataRecords"]:
            raise ValueError(f"Value cell counts do not cover the whole column: {column['key']}")
        column.update({
            "valueCount": len(values), "cellCount": sum(value["cellCount"] for value in values),
            "valueStatusCounts": dict(sorted(Counter(value["status"] for value in values).items())),
            "cellStatusCounts": dict(sorted(sum_counters((value["status"], value["cellCount"]) for value in values).items())),
            "distinctRawValues": len({value["rawValue"] for value in values}),
            "assignedTermCount": len({value["term"]["id"] for value in values if value["term"]}),
            "assignedUnitCount": len({value["unit"]["id"] for value in values if value["unit"]}),
            "termCellCount": sum(value["cellCount"] for value in values if value["term"]),
            "unitCellCount": sum(value["cellCount"] for value in values if value["unit"]),
        })
        # Show frequent values first; the mapping ID is a stable tie-breaker.
        values.sort(key=lambda value: (-value["cellCount"], value["rawValue"], value["mappingId"]))

    projects = []
    for project in sorted({column["project"] for column in columns}):
        selected = [column for column in columns if column["project"] == project]
        project_files = [source for source in files if source["project"] == project]
        projects.append({
            "id": project, "title": titles.get(project, ""),
            "fileCount": len(project_files), "columnCount": len(selected),
            "ontologyColumnCount": sum(column["columnTerm"] is not None for column in selected),
            "reviewColumnCount": sum(column["status"] == "needs_review" for column in selected),
            "cellCount": sum(column["cellCount"] for column in selected),
            "valueCount": sum(column["valueCount"] for column in selected),
            "termCount": len({column["columnTerm"]["id"] for column in selected if column["columnTerm"]}),
            "sourcePaths": sorted(source["sourcePath"] for source in project_files),
        })

    groups = []
    for group_key in sorted({column["groupKey"] for column in columns}):
        selected = [column for column in columns if column["groupKey"] == group_key]
        concept = next((column["columnTerm"] for column in selected if column["columnTerm"]), None)
        groups.append({
            "key": group_key, "label": concept["label"] if concept else selected[0]["standardName"].replace("_", " "),
            "termId": concept["id"] if concept else "",
            "standardNames": sorted({column["standardName"] for column in selected}),
            "projectIds": sorted({column["project"] for column in selected}),
            "columnKeys": [column["key"] for column in selected],
        })

    index = {
        "schemaVersion": 1,
        "meta": {
            "sourceFiles": source_files, "datasetSha256": dataset_sha,
            "counts": {
                "projects": len(projects), "files": len(files), "columns": len(columns),
                "rawHeaders": len({column["sourceColumn"] for column in columns}),
                "standardNames": len({column["standardName"] for column in columns}),
                "ontologyColumns": sum(column["columnTerm"] is not None for column in columns),
                "ontologyTerms": len({column["columnTerm"]["id"] for column in columns if column["columnTerm"]}),
                "valueMappings": len(values_raw), "cells": sum(value["cellCount"] for values in values_by_column.values() for value in values),
                "sourceRecords": sum(source["dataRecords"] for source in files),
            },
            "ontologyColumnStatusCounts": dict(sorted(Counter(column["status"] for column in columns if column["columnTerm"]).items())),
            "notes": [
                "Column ontology categories and cell-value ontology terms describe different things.",
                "A column may retain a known category or explicit unit while its detailed meaning still needs review.",
                "Value mappings are context scoped: identical raw strings can have different mappings or sample kinds.",
                "Cell counts are weighted by cell_count; source records include the header as record 1.",
            ],
        },
        "projects": projects, "columns": columns, "groups": groups,
        "files": sorted(files, key=lambda source: (source["project"], source["sourcePath"])),
    }
    shards = {project["id"]: {"project": project["id"], "columns": {
        column["key"]: values_by_column[column["key"]] for column in columns if column["project"] == project["id"]}}
        for project in projects}
    return index, shards


def sum_counters(pairs):
    result = Counter()
    for key, count in pairs:
        result[key] += count
    return result


def write_snapshot(index, shards, output=DEFAULT_OUTPUT, column_path=None):
    output = Path(output)
    values_dir = output / "values"
    values_dir.mkdir(parents=True, exist_ok=True)
    (output / "explorer.json").write_text(stable_json(index) + "\n", encoding="utf-8")
    expected = {project + ".json" for project in shards}
    for stale in values_dir.glob("*.json"):
        if stale.name not in expected:
            stale.unlink()
    for project, data in shards.items():
        (values_dir / (project + ".json")).write_text(stable_json(data) + "\n", encoding="utf-8")
    if column_path:
        shutil.copyfile(column_path, output / "column_mappings.tsv")
    return {"indexBytes": (output / "explorer.json").stat().st_size,
            "valueBytes": sum(path.stat().st_size for path in values_dir.glob("*.json"))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports-root", type=Path, default=DEFAULT_REPORTS)
    parser.add_argument("--projects-file", type=Path, default=DEFAULT_TITLES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    index, shards = build_snapshot(args.reports_root, args.projects_file)
    sizes = write_snapshot(index, shards, args.output_dir, args.reports_root / "column_mappings.tsv")
    print(json.dumps({"validation": "passed", **index["meta"]["counts"], **sizes,
                      "datasetSha256": index["meta"]["datasetSha256"]}, indent=2))


if __name__ == "__main__":
    main()
