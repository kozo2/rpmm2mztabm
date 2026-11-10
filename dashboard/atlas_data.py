"""Pure data access and query helpers for the RPMM Streamlit dashboard.

Source identities and scientific annotations are retained as recorded. Ontology
categories on columns are distinct from assigned terms on cell values.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import importlib.util
import io
import json
from pathlib import Path
from typing import Iterable

import pandas as pd


DASHBOARD = Path(__file__).resolve().parent
REPOSITORY = DASHBOARD.parent
DEFAULT_REPORTS = REPOSITORY / "reports" / "phenodata"
DEFAULT_SNAPSHOT = DASHBOARD / "snapshot"


@dataclass
class AtlasData:
    index: dict
    shards: dict
    columns: pd.DataFrame
    projects: pd.DataFrame
    values: pd.DataFrame
    source_mode: str
    source_label: str


def _generator():
    path = DASHBOARD / "scripts" / "build-data.py"
    spec = importlib.util.spec_from_file_location("rpmm_dashboard_builder", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _term_field(term, name):
    return term.get(name, "") if isinstance(term, dict) else ""


def _check_snapshot(index, shards):
    if index.get("schemaVersion") != 1:
        raise ValueError("Unsupported dashboard snapshot schema")
    columns = index["columns"]
    keys = [column["key"] for column in columns]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate source-column identity in dashboard data")
    project_ids = {project["id"] for project in index["projects"]}
    if set(shards) != project_ids:
        raise ValueError("Snapshot value projects disagree with the column index")
    mapping_ids = set()
    cells = mappings = 0
    for project in sorted(project_ids):
        shard = shards[project]
        expected = {column["key"] for column in columns if column["project"] == project}
        if shard.get("project") != project or set(shard["columns"]) != expected:
            raise ValueError(f"Values contain missing or foreign source columns: {project}")
        for key, values in shard["columns"].items():
            column = next(column for column in columns if column["key"] == key)
            if sum(value["cellCount"] for value in values) != column["cellCount"]:
                raise ValueError(f"Value counts disagree with source-column counts: {key}")
            for value in values:
                if value["mappingId"] in mapping_ids:
                    raise ValueError("Duplicate source value-mapping identity")
                mapping_ids.add(value["mappingId"])
                cells += value["cellCount"]
                mappings += 1
    counts = index["meta"]["counts"]
    if (len(project_ids) != counts["projects"] or len(columns) != counts["columns"]
            or cells != counts["cells"] or mappings != counts["valueMappings"]):
        raise ValueError("Dashboard snapshot census is inconsistent")


def load_data(reports_root=None, snapshot_dir=None, prefer_reports=True, projects_file=None):
    """Read current report TSVs, or a checked-in snapshot when reports are absent.

    A partial or malformed live report set fails explicitly. It is never quietly
    replaced by older data. `prefer_reports=False` explicitly selects a snapshot.
    """
    reports = Path(reports_root) if reports_root is not None else DEFAULT_REPORTS
    snapshot = Path(snapshot_dir) if snapshot_dir is not None else DEFAULT_SNAPSHOT
    required = [reports / name for name in ("column_mappings.tsv", "value_mappings.tsv", "inventory.tsv")]
    if prefer_reports and any(path.exists() for path in required):
        if not all(path.is_file() for path in required):
            raise FileNotFoundError("The repository report set is incomplete: " + str(reports))
        generator = _generator()
        index, shards = generator.build_snapshot(reports, projects_file or generator.DEFAULT_TITLES)
        mode, label = "reports", "reports/phenodata"
    else:
        index = json.loads((snapshot / "explorer.json").read_text(encoding="utf-8"))
        shards = {project["id"]: json.loads((snapshot / "values" / (project["id"] + ".json")).read_text(encoding="utf-8"))
                  for project in index["projects"]}
        mode, label = "snapshot", "Committed dashboard snapshot"
    _check_snapshot(index, shards)

    titles = {project["id"]: project["title"] for project in index["projects"]}
    column_records = []
    for record in index["columns"]:
        column_records.append({**record, "projectTitle": titles.get(record["project"], ""),
                               "columnTermId": _term_field(record["columnTerm"], "id"),
                               "columnTermLabel": _term_field(record["columnTerm"], "label"),
                               "columnTermIri": _term_field(record["columnTerm"], "iri"),
                               "unitId": _term_field(record["unit"], "id"),
                               "unitLabel": _term_field(record["unit"], "label")})
    columns = pd.DataFrame(column_records)
    column_by_key = {column["key"]: column for column in column_records}
    value_records = []
    for project in sorted(shards):
        for key, records in shards[project]["columns"].items():
            column = column_by_key[key]
            for record in records:
                value_records.append({**record, "columnKey": key, "project": project,
                                      "sourcePath": column["sourcePath"], "sourceColumn": column["sourceColumn"],
                                      "standardName": column["standardName"], "semanticRole": column["semanticRole"],
                                      "columnTermId": column["columnTermId"],
                                      "termId": _term_field(record["term"], "id"),
                                      "termLabel": _term_field(record["term"], "label"),
                                      "termIri": _term_field(record["term"], "iri"),
                                      "unitId": _term_field(record["unit"], "id"),
                                      "unitLabel": _term_field(record["unit"], "label")})
    values = pd.DataFrame(value_records)
    return AtlasData(index=index, shards=shards, columns=columns,
                     projects=pd.DataFrame(index["projects"]), values=values,
                     source_mode=mode, source_label=label)


def _selection(values):
    if values is None:
        return None
    return [values] if isinstance(values, str) else list(values)


def filter_columns(columns, project_ids=None, query="", ontology_ids=None,
                   standard_names=None, semantic_roles=None, statuses=None, mapped_only=False):
    """Filter source-column instances; None ignores a selector, [] selects none."""
    mask = pd.Series(True, index=columns.index)
    for field, selected in (("project", project_ids), ("columnTermId", ontology_ids),
                            ("standardName", standard_names), ("semanticRole", semantic_roles),
                            ("status", statuses)):
        selected = _selection(selected)
        if selected is not None:
            mask &= columns[field].isin(selected)
    if mapped_only:
        mask &= columns["columnTermId"].ne("")
    if query and query.strip():
        searchable = [field for field in ("project", "projectTitle", "sourcePath", "sourceColumn", "standardName",
                                          "semanticRole", "columnTermId", "columnTermLabel", "unitId", "unitLabel")
                      if field in columns]
        text = columns[searchable].fillna("").astype(str).agg(" ".join, axis=1)
        mask &= text.str.contains(query.strip(), case=False, regex=False)
    return columns.loc[mask].copy()


def get_values(data, column_keys=None):
    """Return exact contextual mappings for selected source-column identities."""
    keys = _selection(column_keys)
    return (data.values if keys is None else data.values.loc[data.values["columnKey"].isin(keys)]).copy()


VALUE_GROUP_FIELDS = ("rawValue", "normalizedValue", "numericValue", "termId", "unitId", "status", "rowKind", "referenceEvent")


def _freeze(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) if isinstance(value, (dict, list)) else value


def aggregate_values(values, group_by=None):
    """Aggregate equivalent mappings with cell weights and full source identity.

    Default grouping preserves raw/normalized/numeric values, terms, units,
    status, sample kind, time reference and candidate/assertion differences.
    Different source columns remain separate even when their strings coincide.
    Original mapping IDs and one-based source records remain available.
    """
    fields = list(VALUE_GROUP_FIELDS if group_by is None else group_by)
    for field in fields:
        if field not in values.columns:
            raise KeyError("Unknown value grouping field: " + field)
    if "columnKey" in values and "columnKey" not in fields:
        fields.insert(0, "columnKey")
    for field in ("candidateTermId", "candidateUnitId", "candidateTerms", "assertions"):
        if field in values and field not in fields:
            fields.append(field)
    if values.empty:
        return pd.DataFrame(columns=list(dict.fromkeys([*values.columns, "mappingCount", "mappingIds"])))
    buckets = {}
    for row in values.to_dict("records"):
        key = tuple(_freeze(row.get(field)) for field in fields)
        if key not in buckets:
            record = {name: value for name, value in row.items() if name != "mappingId"}
            record.update(cellCount=0, mappingCount=0, mappingIds=[], sourceRecords=[], evidenceVariants=[], reasonVariants=[])
            buckets[key] = record
        record = buckets[key]
        record["cellCount"] += int(row["cellCount"])
        record["mappingCount"] += 1
        if row.get("mappingId"):
            record["mappingIds"].append(row["mappingId"])
        record["sourceRecords"].extend(row.get("sourceRecords", []))
        record["evidenceVariants"].append(row.get("evidence", ""))
        record["reasonVariants"].append(row.get("reason", ""))
    records = list(buckets.values())
    for record in records:
        for field in ("mappingIds", "sourceRecords", "evidenceVariants", "reasonVariants"):
            record[field] = sorted(set(record[field]))
    records.sort(key=lambda row: (-row["cellCount"], str(row.get("rawValue", "")), str(row.get("columnKey", ""))))
    return pd.DataFrame(records)


def build_matrix(columns, group_by="ontology", project_ids=None):
    """Count source-column instances by project and category or standard field.

    Ontology matrix columns are actual category IDs, not labels. Projects can
    explicitly be reindexed to include zero rows; cells count column instances,
    not samples, assignments, or booleans.
    """
    if group_by not in {"ontology", "standard"}:
        raise ValueError("group_by must be ontology or standard")
    selected_projects = _selection(project_ids)
    selected_projects = sorted(columns["project"].unique()) if selected_projects is None else list(dict.fromkeys(selected_projects))
    field = "columnTermId" if group_by == "ontology" else "standardName"
    chosen = columns.loc[columns["project"].isin(selected_projects)]
    chosen = chosen.loc[chosen[field].ne("")]
    if chosen.empty:
        return pd.DataFrame(index=pd.Index(selected_projects, name="project"), dtype="int64")
    matrix = chosen.groupby(["project", field], sort=True).size().unstack(fill_value=0)
    matrix = matrix.reindex(index=selected_projects, fill_value=0)
    matrix.index.name = "project"
    matrix.columns.name = field
    return matrix.astype("int64")


def to_csv(frame):
    """Export UTF-8 CSV without pandas NA/string coercion or lossy numbers."""
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(frame.columns)
    for row in frame.itertuples(index=False, name=None):
        rendered = []
        for value in row:
            if isinstance(value, (dict, list)):
                rendered.append(json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")))
            elif value is None or (not isinstance(value, str) and pd.isna(value)):
                rendered.append("")
            else:
                rendered.append(value)
        writer.writerow(rendered)
    return output.getvalue().encode("utf-8")
