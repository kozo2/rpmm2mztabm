"""Regression cases for the dashboard's scientific joins, counts and exports."""

import csv
import io
import json
from pathlib import Path
import sys
import unittest

import pandas as pd

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))

from atlas_data import aggregate_values, build_matrix, filter_columns, get_values, load_data, to_csv


class AtlasQueriesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_data()

    def columns(self, project, header=None):
        selected = filter_columns(self.data.columns, project_ids=[project])
        return selected if header is None else selected[selected["sourceColumn"] == header]

    def test_authoritative_inventory_and_column_identity(self):
        data = self.data
        self.assertEqual(len(data.projects), 31)
        self.assertEqual(len(data.columns), 391)
        self.assertEqual(data.columns["sourcePath"].nunique(), 39)
        self.assertEqual(data.columns["sourceColumn"].nunique(), 165)
        self.assertEqual(data.columns["key"].nunique(), 391)
        for row in data.columns.to_dict("records"):
            self.assertEqual(row["key"], "|".join((row["project"], row["sourcePath"], str(row["columnIndex"]))))

    def test_snapshot_fallback_preserves_the_same_scientific_data(self):
        snapshot = load_data(prefer_reports=False)
        self.assertEqual(snapshot.index["meta"]["datasetSha256"], self.data.index["meta"]["datasetSha256"])
        self.assertEqual(set(snapshot.columns["key"]), set(self.data.columns["key"]))
        self.assertEqual(set(snapshot.values["mappingId"]), set(self.data.values["mappingId"]))

    def test_column_categories_are_distinct_from_annotation_status(self):
        mapped = filter_columns(self.data.columns, mapped_only=True)
        self.assertEqual(len(mapped), 143)
        self.assertEqual(mapped["status"].value_counts().to_dict(),
                         {"accepted": 100, "local": 23, "needs_review": 20})
        genotype = self.columns("RPMM0068", "genotype")
        self.assertEqual(genotype.iloc[0]["status"], "local")
        values = get_values(self.data, genotype["key"].tolist())
        self.assertEqual(set(values["status"]), {"local"})
        self.assertFalse(values["termId"].astype(bool).any())

    def test_filters_are_case_insensitive_and_explicit_empty_means_empty(self):
        anatomy = filter_columns(self.data.columns, query="efo:0000635")
        self.assertEqual(len(anatomy), 21)
        self.assertTrue(filter_columns(self.data.columns, project_ids=[]).empty)
        self.assertTrue(filter_columns(self.data.columns, query="no_such_ontology_or_column_987654").empty)
        selected = filter_columns(self.data.columns, project_ids=["RPMM0068"], statuses=["needs_review"])
        self.assertEqual(len(selected), 4)

    def test_same_header_in_two_files_is_never_the_same_column(self):
        treatment = self.columns("RPMM0027", "treatment")
        self.assertEqual(len(treatment), 2)
        self.assertEqual(treatment["sourcePath"].nunique(), 2)
        values = get_values(self.data, treatment["key"].tolist())
        acetate = values[values["sourcePath"].str.contains("AcetateTreatment")]
        drought = values[values["sourcePath"].str.contains("DroughtTreatment")]
        self.assertEqual(acetate.groupby("rawValue")["cellCount"].sum().to_dict(), {"AcOH": 5, "Water": 5})
        self.assertEqual(drought.groupby("rawValue")["cellCount"].sum().to_dict(), {"dry": 30, "water": 39})

    def test_repeated_context_mappings_use_weighted_occurrence_counts(self):
        concentration = self.columns("RPMM0068", "conc")
        values = get_values(self.data, concentration["key"].tolist())
        self.assertGreater(len(values), 2)
        aggregated = aggregate_values(values)
        self.assertEqual(aggregated.groupby("rawValue")["cellCount"].sum().to_dict(), {"0": 24, "200": 24})
        self.assertEqual(int(aggregated["cellCount"].sum()), 48)
        self.assertEqual(int(aggregated["mappingCount"].sum()), len(values))

    def test_null_genotype_is_a_literal_under_review_not_missing(self):
        columns = self.columns("RPMM0063", "Genotype2")
        values = get_values(self.data, columns["key"].tolist())
        nulls = values[values["rawValue"] == "null"]
        self.assertEqual(int(nulls["cellCount"].sum()), 42)
        self.assertEqual(set(nulls["status"]), {"needs_review"})
        self.assertEqual(set(nulls["normalizedValue"]), {"null"})

    def test_standard_reference_and_experimental_rows_stay_distinct(self):
        columns = self.columns("RPMM0017", "tissue")
        values = get_values(self.data, columns["key"].tolist())
        leaves = aggregate_values(values[values["rawValue"] == "leaf"])
        self.assertEqual(leaves.groupby("rowKind")["cellCount"].sum().to_dict(),
                         {"experimental_record": 18, "reference_material": 9})
        missing = values[values["rawValue"] == "NA"]
        self.assertEqual(int(missing["cellCount"].sum()), 45)
        self.assertEqual(set(missing["status"]), {"missing"})
        self.assertEqual(set(missing["rowKind"]), {"standard_material"})

    def test_matrix_counts_source_instances_and_preserves_zero_rows(self):
        mapped = filter_columns(self.data.columns, mapped_only=True)
        matrix = build_matrix(mapped, group_by="ontology", project_ids=["RPMM0068", "RPMM0067", "RPMM9999"])
        self.assertEqual(int(matrix.loc["RPMM0068", "PATO:0000125"]), 5)
        self.assertEqual(int(matrix.loc["RPMM0067", "EFO:0000635"]), 3)
        self.assertEqual(int(matrix.loc["RPMM9999"].sum()), 0)
        standards = build_matrix(mapped, group_by="standard", project_ids=["RPMM0068"])
        self.assertEqual(int(standards.loc["RPMM0068", "pooled_plant_mass"]), 2)
        self.assertEqual(int(standards.loc["RPMM0068", "mass_per_plant"]), 2)

    def test_candidate_taxon_ids_are_not_assigned_terms(self):
        columns = self.columns("RPMM0037", "sepices")
        values = get_values(self.data, columns["key"].tolist())
        bacteria = values[values["rawValue"].str.startswith("Bacillus")]
        self.assertFalse(bacteria.empty)
        self.assertEqual(set(bacteria["status"]), {"needs_review"})
        self.assertFalse(bacteria["termId"].astype(bool).any())
        self.assertTrue(any(component.get("candidate_term_id") == "NCBITaxon:1402"
                            for assertions in bacteria["assertions"] for component in assertions))

    def test_complete_weighted_status_totals(self):
        values = get_values(self.data)
        self.assertEqual(len(values), 11784)
        self.assertEqual(int(values["cellCount"].sum()), 26997)
        self.assertEqual(values.groupby("status")["cellCount"].sum().to_dict(),
                         {"accepted": 7567, "local": 14671, "missing": 1051,
                          "needs_review": 3627, "not_applicable": 81})

    def test_csv_roundtrip_preserves_literals_precision_and_structured_candidates(self):
        frame = pd.DataFrame({"rawValue": ["null", "", 'line 1\nline 2,"quoted"'],
                              "numericValue": ["0.1234567890123456789012345", "", "-2"],
                              "assertions": [[{"candidate_term_id": "NCBITaxon:1402"}], [], []]})
        exported = to_csv(frame)
        self.assertIsInstance(exported, bytes)
        rows = list(csv.DictReader(io.StringIO(exported.decode("utf-8-sig"))))
        self.assertEqual([row["rawValue"] for row in rows], frame["rawValue"].tolist())
        self.assertEqual(rows[0]["numericValue"], "0.1234567890123456789012345")
        self.assertEqual(json.loads(rows[0]["assertions"]), frame.iloc[0]["assertions"])


if __name__ == "__main__":
    unittest.main()
