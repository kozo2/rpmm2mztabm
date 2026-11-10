"""Regression checks for provenance, fail-closed IDs, and biological semantics.

Run from the repository root:
  PYTHONDONTWRITEBYTECODE=1 python3 -m unittest investigation.phenodata_ontology.test_annotation

No export is generated and no network access is required.
"""

import hashlib
import csv
import json
import os
from pathlib import Path
import tempfile
import unittest

from . import catalog
from .columns import classify_column
from .pipeline import (
    DEFAULT_DATA,
    add_candidates,
    classify_row,
    detect_conflicts,
    load_inputs,
    validate_annotation,
    validate_column,
)
from .values import annotate_value
from .validate_export import expected_role_coverage, validate_project_exports


def annotated(project, header, value, row=None, path=""):
    """Exercise the same schema/value/catalog boundary as a normal export."""
    spec = validate_column(classify_column(project, path, header))
    result = annotate_value(project, path, header, value,
                            row if row is not None else {header: value}, spec)
    return validate_annotation(result), spec


class SourceInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = load_inputs(DEFAULT_DATA)

    def test_all_source_records_and_cells_are_retained(self):
        self.assertEqual(len(self.inputs), 39)
        self.assertEqual(len({item["project"] for item in self.inputs}), 31)
        self.assertEqual(sum(len(item["headers"]) for item in self.inputs), 391)
        self.assertEqual(sum(len(item["records"]) for item in self.inputs), 2606)
        self.assertEqual(sum(len(record["cells"]) for item in self.inputs
                             for record in item["records"]), 26997)
        nonempty = sum(any(value.strip() for value in record["cells"])
                       for item in self.inputs for record in item["records"])
        self.assertEqual(nonempty, 2602)
        for item in self.inputs:
            with self.subTest(path=item["relative_path"]):
                self.assertEqual(item["sha256"], hashlib.sha256(item["path"].read_bytes()).hexdigest())
                self.assertEqual([record["source_record"] for record in item["records"]],
                                 list(range(2, 2 + len(item["records"]))))
                self.assertTrue(all(len(record["cells"]) == len(item["headers"])
                                    for record in item["records"]))

    def test_original_header_spelling_and_whitespace_survive(self):
        item = next(item for item in self.inputs
                    if "RPMM0012/SupplementalData1/T3/" in item["relative_path"])
        self.assertIn("Phenotype information ", item["headers"])
        self.assertNotIn("Phenotype information", item["headers"])
        # Typographic fixes belong in mappings, never in the input identity.
        spellings = {header for item in self.inputs for header in item["headers"]}
        self.assertIn("Speacies", spellings)
        self.assertIn("sepices", spellings)

    def test_quoted_multiline_header_and_value_have_physical_line_ranges(self):
        raw = (b'gcid,"note\r\nheader"\r\n'
               b'sample1,"line 1\r\nline 2, with comma"\r\n'
               b'sample2,"  preserve spaces  "\r\n')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / "RPMM9999" / "nested" / "Phenodata.csv"
            path.parent.mkdir(parents=True)
            path.write_bytes(raw)
            item = load_inputs(root)[0]
        self.assertEqual(item["headers"], ["gcid", "note\r\nheader"])
        self.assertEqual(item["records"][0]["cells"], ["sample1", "line 1\r\nline 2, with comma"])
        self.assertEqual(item["records"][1]["cells"], ["sample2", "  preserve spaces  "])
        self.assertEqual([(record["source_record"], record["source_line_start"], record["source_line_end"])
                          for record in item["records"]], [(2, 3, 4), (3, 5, 5)])
        self.assertEqual(item["sha256"], hashlib.sha256(raw).hexdigest())

    def test_ambiguous_or_malformed_csv_is_rejected(self):
        for raw in (b'gcid,gcid\nfirst,second\n', b'gcid,note\nfirst,second,extra\n'):
            with self.subTest(raw=raw), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                path = root / "RPMM9999" / "phenodata.csv"
                path.parent.mkdir()
                path.write_bytes(raw)
                with self.assertRaises(ValueError):
                    load_inputs(root)

    def test_project_data_layout_preserves_logical_paths_and_ignores_derivatives(self):
        raw = b'gcid,tissue\nsample1,leaf\n'
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            old = root / "legacy" / "RPMM9999" / "AnalysedDataset" / "phenodata.csv"
            new = root / "projects" / "RPMM9999" / "data" / "AnalysedDataset" / "phenodata.csv"
            derivative = root / "projects" / "RPMM9999" / "phenodata" / "phenodata.csv"
            for path in (old, new, derivative):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(raw)
            legacy = load_inputs(root / "legacy")
            relocated = load_inputs(root / "projects")
            self.assertEqual(len(relocated), 1)
            for key in ("relative_path", "project", "sha256", "md5", "headers", "records"):
                self.assertEqual(relocated[0][key], legacy[0][key])
            self.assertEqual(relocated[0]["relative_path"], "RPMM9999/AnalysedDataset/phenodata.csv")

    def test_default_source_root_is_independent_of_working_directory(self):
        previous = Path.cwd()
        with tempfile.TemporaryDirectory() as folder:
            try:
                os.chdir(folder)
                self.assertEqual(len(load_inputs(DEFAULT_DATA)), 39)
            finally:
                os.chdir(previous)


class ProjectExportTests(unittest.TestCase):
    @staticmethod
    def _write(path, rows):
        fields = sorted({key for row in rows for key in row}) or ["project"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def _fixture(self, folder):
        root = Path(folder)
        output = root / "reports" / "phenodata"
        projects = root / "projects"
        ids = ["RPMM0003", "RPMM0004"]
        tables = {
            "inventory": [{"project": p, "source_path": p + "/AnalysedDataset/phenodata.csv"} for p in ids],
            "cell_annotations": [{"project": p, "semantic_role": "organism", "status": "accepted",
                                  "term_id": "NCBITaxon:3702", "unit_id": ""} for p in ids],
            "value_mappings": [{"project": p, "semantic_role": "organism", "status": "accepted",
                                "term_id": "NCBITaxon:3702", "unit_id": ""} for p in ids],
            "column_mappings": [{"project": p, "source_column": "species"} for p in ids],
            "row_annotations": [{"project": p, "row_kind": "experimental_record"} for p in ids],
            "curation_queue": [{"project": "", "issue_code": "external_source"},
                               *[{"project": p, "issue_code": "review"} for p in ids]],
        }
        report = {"rule_file_sha256": {"rules.py": "rules-hash"},
                  "catalog_file_sha256": {"terms.json": "catalog-hash"},
                  "evidence_file_sha256": {"manifest.csv": "evidence-hash"},
                  "project_exports": {p: {"directory": os.path.relpath(projects / p / "phenodata", output)}
                                      for p in ids}}
        self._write(output / "coverage_by_project.tsv", [{"project": p, "cells": "1"} for p in ids])
        (output / "qa_report.json").write_text(json.dumps(report))
        for project in ids:
            directory = projects / project / "phenodata"
            selected = {name: [row for row in rows if row["project"] == project]
                        for name, rows in tables.items()}
            for name, rows in selected.items():
                self._write(directory / (name + ".tsv"), rows)
            self._write(directory / "coverage_by_project.tsv", [{"project": project, "cells": "1"}])
            self._write(directory / "coverage_by_role.tsv",
                        expected_role_coverage(selected["cell_annotations"], selected["value_mappings"]))
            provenance = {key: report[key] for key in (
                "rule_file_sha256", "catalog_file_sha256", "evidence_file_sha256")}
            provenance.update(project=project,
                              aggregate_report=os.path.relpath(output / "qa_report.json", directory),
                              files=1, source_cells=1, source_records=1,
                              source_column_instances=1, unique_value_mappings=1)
            (directory / "provenance.json").write_text(json.dumps(provenance))
        return output, projects, report, tables

    def test_project_shards_equal_aggregate_and_omit_global_only_rows(self):
        with tempfile.TemporaryDirectory() as folder:
            output, projects, report, tables = self._fixture(folder)
            self.assertEqual(validate_project_exports(output, report, tables), 2)
            self.assertEqual(validate_project_exports(output, report, tables, projects), 2)

    def test_project_shard_cannot_contain_another_project(self):
        with tempfile.TemporaryDirectory() as folder:
            output, projects, report, tables = self._fixture(folder)
            self._write(projects / "RPMM0003" / "phenodata" / "cell_annotations.tsv",
                        tables["cell_annotations"])
            with self.assertRaisesRegex(ValueError, "Another project's rows"):
                validate_project_exports(output, report, tables)

    def test_project_shard_cannot_silently_omit_rows_or_change_provenance(self):
        for defect in ("missing_row", "changed_hash", "broken_reference"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as folder:
                output, projects, report, tables = self._fixture(folder)
                directory = projects / "RPMM0003" / "phenodata"
                if defect == "missing_row":
                    self._write(directory / "value_mappings.tsv", [])
                else:
                    path = directory / "provenance.json"
                    provenance = json.loads(path.read_text())
                    if defect == "changed_hash":
                        provenance["catalog_file_sha256"] = {"terms.json": "different"}
                    else:
                        provenance["aggregate_report"] = "missing.json"
                    path.write_text(json.dumps(provenance))
                with self.assertRaises(ValueError):
                    validate_project_exports(output, report, tables)


class CatalogBoundaryTests(unittest.TestCase):
    def test_invalid_later_component_downgrades_earlier_verified_components(self):
        result = validate_annotation({
            "term_id": "PO:0009005", "status": "accepted",
            "assertions": [{"term_id": "PO:0025034"}, {"term_id": "PO:9999999"}],
        })
        self.assertEqual(result["status"], "needs_review")
        self.assertFalse(result["term_id"])
        for assertion in result["assertions"]:
            self.assertFalse(assertion.get("term_id"))
            self.assertTrue(assertion.get("candidate_term_id"))

    def test_obsolete_and_unknown_ids_are_candidates_not_assignments(self):
        for term_id in ("OBI:0100026", "PO:9999999"):
            with self.subTest(term_id=term_id):
                result = validate_annotation({"term_id": term_id, "status": "accepted"})
                self.assertEqual(result["status"], "needs_review")
                self.assertEqual(result["term_id"], "")
                self.assertEqual(result["candidate_term_id"], term_id)
        self.assertFalse(catalog.is_current("OBI:0100026"))
        self.assertEqual(catalog.get("OBI:0100026")["replaced_by"], ["COB:0000022"])
        self.assertTrue(catalog.is_current("COB:0000022"))

    def test_unknown_unit_is_also_retained_as_candidate(self):
        result = validate_annotation({"numeric_value": "5", "unit_id": "UO:9999999", "status": "accepted"})
        self.assertEqual(result["numeric_value"], "5")
        self.assertEqual(result["unit_id"], "")
        self.assertEqual(result["candidate_unit_id"], "UO:9999999")
        self.assertEqual(result["status"], "needs_review")

    def test_review_cannot_leak_an_assigned_id_from_components(self):
        result = validate_annotation({
            "term_id": "PO:0009005", "unit_id": "UO:0000022", "status": "needs_review",
            "assertions": [{"role": "candidate_structure", "term_id": "PO:0025034"}],
        })
        self.assertEqual(result["term_id"], "")
        self.assertEqual(result["unit_id"], "")
        self.assertEqual(result["candidate_term_id"], "PO:0009005")
        self.assertEqual(result["candidate_unit_id"], "UO:0000022")
        self.assertEqual(result["assertions"][0]["term_id"], "")
        self.assertEqual(result["assertions"][0]["candidate_term_id"], "PO:0025034")

    def test_current_assignment_carries_native_label_version_and_source(self):
        result = validate_annotation({"term_id": "PO:0025034", "status": "accepted"})
        self.assertEqual(result["term_label"], "leaf")
        self.assertTrue(result["term_iri"].endswith("PO_0025034"))
        self.assertTrue(result["term_source_id"])
        self.assertTrue(result["term_ontology_version"])

    def test_critical_ids_match_intended_species_and_anatomy(self):
        expected = {
            "NCBITaxon:35875": "Allium fistulosum",
            "NCBITaxon:130426": "Allium chinense",
            "NCBITaxon:4683": "Allium tuberosum",
            "NCBITaxon:210367": "Thesium chinense",
            "PO:0009072": "plant ovary",
            "PO:0009084": "pericarp",
            "PATO:0000125": "mass",
            "UO:0000022": "milligram",
        }
        for term_id, label in expected.items():
            with self.subTest(term_id=term_id):
                self.assertTrue(catalog.is_current(term_id))
                self.assertEqual(catalog.get(term_id)["label"], label)

    def test_native_candidate_is_retained_without_being_autoapproved(self):
        result = add_candidates({"status": "needs_review", "normalized_value": "root"},
                                {"semantic_role": "organism_part"}, "root")
        self.assertEqual(result["status"], "needs_review")
        self.assertFalse(result.get("term_id"))
        self.assertIn("PO:0009005", {term["term_id"] for term in result["candidate_terms"]})


class BiologicalSemanticsTests(unittest.TestCase):
    def test_standards_and_references_are_not_experimental_plant_rows(self):
        for row in ({"sample": "standard"}, {"sample": "Standard mixture"},
                    {"sample": "Standards mixture"}, {"class1": "STDs"},
                    {"type": "SD"}, {"type": "alkanemix"}):
            with self.subTest(row=row):
                self.assertEqual(classify_row(row)[0], "standard_material")
        self.assertEqual(classify_row({"fruitvsleaf": "Freeze Dry Col-0"})[0], "reference_material")
        self.assertEqual(classify_row({"gcid": "first", "species": "tomato"})[0], "experimental_record")
        self.assertEqual(classify_row({"gcid": "", "species": ""})[0], "empty")

    def test_standard_markers_do_not_become_taxa_or_tissues(self):
        for header in ("species", "tissue"):
            for value in ("STDs", "standard", "Standards mixture"):
                with self.subTest(header=header, value=value):
                    result, _ = annotated("RPMM0003", header, value)
                    self.assertFalse(result["term_id"])
                    self.assertEqual(result["status"], "not_applicable")
            result, _ = annotated("RPMM0003", header, "Freeze Dry Col-0")
            self.assertFalse(result["term_id"])
            # A biological reference still has biological properties; this
            # marker does not state their values or make them inapplicable.
            self.assertEqual(result["status"], "needs_review")

    def test_null_genotype_is_preserved_and_distinct_from_missing(self):
        result, _ = annotated("RPMM0063", "Genotype2", "null")
        self.assertEqual(result["value_kind"], "identifier")
        self.assertEqual(result["normalized_value"], "null")
        self.assertEqual(result["status"], "needs_review")
        for missing in ("", "NA", "nodata", "noinfo"):
            with self.subTest(missing=missing):
                result, _ = annotated("RPMM0063", "Genotype2", missing)
                self.assertEqual(result["status"], "missing")

    def test_col0_is_local_line_identity_not_taxonomy_or_stock_accession(self):
        for spelling in ("Col-0", "col0", "Col0"):
            with self.subTest(spelling=spelling):
                result, _ = annotated("RPMM0019", "genotype", spelling)
                self.assertEqual(result["normalized_value"], "Col-0")
                self.assertEqual(result["status"], "local")
                self.assertEqual(result["value_kind"], "identifier")
                self.assertFalse(result["term_id"])

    def test_treatment_chemical_has_agent_scope(self):
        result, _ = annotated("RPMM0027", "treatment", "AcOH")
        self.assertEqual(result["term_id"], "CHEBI:15366")
        self.assertEqual(result["term_label"], "acetic acid")
        self.assertTrue(any(assertion.get("role") == "annotation_scope"
                            and assertion.get("value") == "chemical_agent_only"
                            for assertion in result["assertions"]))

    def test_flowering_before_after_offsets_keep_sign_and_reference(self):
        for value, expected in (("2dbf", "-2"), ("0daf", "0"), ("4daf", "4")):
            with self.subTest(value=value):
                result, spec = annotated("RPMM0064", "floweringday", value)
                self.assertEqual(spec["column_term_id"], "EFO:0000724")
                self.assertEqual(result["numeric_value"], expected)
                self.assertEqual(result["unit_id"], "UO:0000033")
                self.assertEqual(result["reference_event"], "flowering")
                self.assertEqual(result["status"], "accepted")

    def test_unknown_quantity_units_are_not_inferred_from_magnitude(self):
        result, _ = annotated("RPMM0010", "incubationtemperature", "60")
        self.assertEqual(result["numeric_value"], "60")
        self.assertFalse(result["unit_id"])
        self.assertEqual(result["status"], "needs_review")
        result, _ = annotated("RPMM0068", "plantweight_gc", "50.9")
        self.assertFalse(result["unit_id"])
        self.assertEqual(result["status"], "needs_review")

    def test_column_label_day_does_not_override_explicit_week_suffix(self):
        result, _ = annotated("RPMM0067", "day_after_cutting", "4w")
        self.assertEqual(result["numeric_value"], "4")
        self.assertEqual(result["unit_id"], "UO:0000034")
        self.assertEqual(result["reference_event"], "cutting")
        self.assertEqual(result["status"], "accepted")

    def test_same_row_genus_and_epithet_determine_allium_species(self):
        for epithet, term_id in (("fistulosum", "NCBITaxon:35875"),
                                 ("chinense", "NCBITaxon:130426"),
                                 ("tuberosum", "NCBITaxon:4683")):
            with self.subTest(epithet=epithet):
                row = {"material": "allium", "speciesname": "spring onion", "binomialname": epithet}
                result, _ = annotated("RPMM0004", "binomialname", epithet, row)
                self.assertEqual(result["term_id"], term_id)
        result, _ = annotated("RPMM0004", "material", "allium", {"material": "allium"},
                              path="RPMM0004/ITEX/AnalysedDataset/phenodata.csv")
        self.assertEqual(result["term_id"], "NCBITaxon:4678")
        result, _ = annotated("RPMM0004", "binomialname", "fistulosum", {"material": "unknown"})
        self.assertFalse(result["term_id"])
        self.assertEqual(result["status"], "needs_review")

    def test_genotype_column_in_fermentation_is_not_a_plant_line(self):
        spec = classify_column("RPMM0037", "RPMM0037/AnalysedDataset/phenodata.csv", "genotype")
        self.assertNotEqual(spec["semantic_role"], "genotype")
        result, _ = annotated("RPMM0037", "sepices", "Bacillus licheniformis 67")
        self.assertEqual(result["status"], "needs_review")
        self.assertFalse(result["term_id"])
        self.assertTrue(any(assertion.get("candidate_term_id") == "NCBITaxon:1402"
                            for assertion in result["assertions"]))

    def test_mixture_ratio_conflict_is_detected_without_choosing_a_side(self):
        item = {"project": "RPMM0017", "relative_path": "RPMM0017/AnalysedDataset/valid02_gc/Phenodata.csv"}
        row = {"mixtureratiofruit": "100", "mixtureratioleaf": "0", "note": "leaf : fruit = 100 : 0"}
        issues = detect_conflicts(item, 2, row)
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0]["issue_code"], "conflicting_mixture_ratio")
        self.assertEqual(set(issues[0]["columns"]), {"mixtureratiofruit", "mixtureratioleaf", "note"})
        consistent = {**row, "note": "leaf : fruit = 0 : 100"}
        self.assertEqual(detect_conflicts(item, 2, consistent), [])


if __name__ == "__main__":
    unittest.main()
