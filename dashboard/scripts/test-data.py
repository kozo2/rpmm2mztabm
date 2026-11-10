#!/usr/bin/env python3
"""Data integrity regressions for the dashboard's scientific source joins."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).with_name("build-data.py")
spec = importlib.util.spec_from_file_location("dashboard_data", SCRIPT)
data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data)


class DashboardDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index, cls.shards = data.build_snapshot()
        cls.source_values = data.read_tsv(data.DEFAULT_REPORTS / "value_mappings.tsv")

    def test_all_source_columns_values_and_cells_are_accounted_for(self):
        counts = self.index["meta"]["counts"]
        self.assertEqual({key: counts[key] for key in ("projects", "files", "columns", "valueMappings", "cells")},
                         {"projects": 31, "files": 39, "columns": 391, "valueMappings": 11784, "cells": 26997})
        columns = self.index["columns"]
        self.assertEqual(len({column["key"] for column in columns}), len(columns))
        self.assertEqual(sum(column["cellCount"] for column in columns), counts["cells"])
        self.assertEqual(sum(project["cellCount"] for project in self.index["projects"]), counts["cells"])

    def test_every_mapping_preserves_source_values_and_candidate_assertions(self):
        values = {value["mappingId"]: value for shard in self.shards.values()
                  for group in shard["columns"].values() for value in group}
        self.assertEqual(len(values), len(self.source_values))
        for raw in self.source_values:
            value = values[raw["mapping_id"]]
            self.assertEqual(value["rawValue"], raw["raw_value"])
            self.assertEqual(value["numericValue"], raw["numeric_value"])
            self.assertEqual(value["normalizedValue"], raw["normalized_value"])
            self.assertEqual(value["cellCount"], int(raw["cell_count"]))
            self.assertEqual(value["sourceRecords"], json.loads(raw["source_records"]))
            self.assertEqual(value["assertions"], json.loads(raw["assertions"] or "[]"))
            self.assertEqual(value["candidateTerms"], json.loads(raw["candidate_terms"] or "[]"))
            self.assertEqual(value["candidateTermId"], raw["candidate_term_id"])
            self.assertEqual(value["candidateUnitId"], raw["candidate_unit_id"])

    def test_file_scope_and_column_index_are_part_of_column_identity(self):
        for column in self.index["columns"]:
            self.assertEqual(column["key"], f"{column['project']}|{column['sourcePath']}|{column['columnIndex']}")
        gcid = [column for column in self.index["columns"]
                if column["project"] == "RPMM0004" and column["sourceColumn"] == "gcid"]
        self.assertEqual(len(gcid), 2)
        self.assertEqual(len({column["key"] for column in gcid}), 2)

    def test_contextual_duplicates_are_weighted_by_cell_count(self):
        column = next(column for column in self.index["columns"]
                      if column["project"] == "RPMM0068" and column["sourceColumn"] == "conc")
        values = self.shards["RPMM0068"]["columns"][column["key"]]
        for raw in ("0", "200"):
            selected = [value for value in values if value["rawValue"] == raw]
            self.assertGreater(len(selected), 1)
            self.assertEqual(sum(value["cellCount"] for value in selected), 24)
        self.assertEqual(column["cellCount"], 48)
        self.assertEqual(column["distinctRawValues"], 2)

    def test_ontology_category_does_not_collapse_distinct_standard_fields(self):
        mass = next(group for group in self.index["groups"] if group["termId"] == "PATO:0000125")
        self.assertGreater(len(mass["standardNames"]), 1)
        self.assertIn("recorded_mass", mass["standardNames"])
        self.assertIn("mass_per_plant", mass["standardNames"])
        self.assertEqual(self.index["meta"]["counts"]["ontologyColumns"], 143)
        self.assertEqual(self.index["meta"]["ontologyColumnStatusCounts"],
                         {"accepted": 100, "local": 23, "needs_review": 20})

    def test_null_genotype_is_preserved_and_review_candidates_are_not_assignments(self):
        nulls = [value for values in self.shards["RPMM0063"]["columns"].values()
                 for value in values if value["rawValue"] == "null"]
        self.assertEqual(sum(value["cellCount"] for value in nulls), 42)
        self.assertTrue(all(value["status"] == "needs_review" for value in nulls))
        for shard in self.shards.values():
            for values in shard["columns"].values():
                for value in values:
                    if value["status"] == "needs_review":
                        self.assertIsNone(value["term"])
                        self.assertIsNone(value["unit"])

    def test_multiline_header_and_trailing_whitespace_are_not_rewritten(self):
        names = {column["sourceColumn"] for column in self.index["columns"]}
        self.assertIn("sample information\n_Genotype_ Phenotype_ Tissue", names)
        self.assertIn("Phenotype information ", names)

    def test_project_titles_are_plain_text_without_losing_entities(self):
        self.assertEqual(data.title_text("A <i>plant</i> &amp; its <sup>2</sup> leaves"), "A plant & its 2 leaves")
        self.assertTrue(all("<i>" not in project["title"] for project in self.index["projects"]))

    def test_generated_snapshot_is_deterministic_and_download_is_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            first, second = Path(temporary) / "first", Path(temporary) / "second"
            for folder in (first, second):
                data.write_snapshot(self.index, self.shards, folder, data.DEFAULT_REPORTS / "column_mappings.tsv")
            for original in first.rglob("*"):
                if original.is_file():
                    self.assertEqual(original.read_bytes(), (second / original.relative_to(first)).read_bytes())
            self.assertEqual((first / "column_mappings.tsv").read_bytes(),
                             (data.DEFAULT_REPORTS / "column_mappings.tsv").read_bytes())
        sources = self.index["meta"]["sourceFiles"]
        for source in sources:
            self.assertEqual(hashlib.sha256((data.REPOSITORY / source["path"]).read_bytes()).hexdigest(), source["sha256"])

    def test_orphan_mapping_is_rejected_instead_of_silently_joined_by_header(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            for name in ("column_mappings.tsv", "inventory.tsv"):
                (folder / name).write_bytes((data.DEFAULT_REPORTS / name).read_bytes())
            raw = dict(self.source_values[0])
            raw["column_index"] = "999"
            with (folder / "value_mappings.tsv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(raw), delimiter="\t")
                writer.writeheader()
                writer.writerow(raw)
            with self.assertRaisesRegex(ValueError, "nonexistent source column"):
                data.build_snapshot(folder)


if __name__ == "__main__":
    unittest.main()
