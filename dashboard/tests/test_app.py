"""Check the focused project -> source column -> recorded values workflow."""

from pathlib import Path
import sys
import unittest

import pandas as pd
from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))

from atlas_data import filter_columns, load_data


class AtlasAppTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_data()

    def app(self):
        app = AppTest.from_file(str(APP / "app.py"), default_timeout=20).run()
        self.healthy(app)
        return app

    def healthy(self, app):
        self.assertFalse(app.exception, [error.message for error in app.exception])

    def key(self, project, header, source_contains=""):
        columns = filter_columns(self.data.columns, project_ids=[project])
        columns = columns[columns["sourceColumn"] == header]
        if source_contains:
            columns = columns[columns["sourcePath"].str.contains(source_contains, regex=False)]
        self.assertEqual(len(columns), 1)
        return columns.iloc[0]["key"]

    def inspect(self, app, project, header, source_contains=""):
        app.selectbox(key="project").set_value(project).run()
        self.healthy(app)
        app.selectbox(key="column").set_value(self.key(project, header, source_contains)).run()
        self.healthy(app)
        return self.value_frame(app)

    def value_frame(self, app):
        matches = [element.value for element in app.dataframe
                   if {"Value", "Occurrences"}.issubset(element.value.columns)]
        return matches[0] if matches else pd.DataFrame(columns=["Value", "Occurrences"])

    def column_frame(self, app):
        return next(element.value for element in app.dataframe
                    if {"Column", "Ontology ID", "Source file"}.issubset(element.value.columns))

    def counts(self, frame):
        return frame.groupby("Value", dropna=False)["Occurrences"].sum().to_dict()

    def test_default_view_is_a_single_project_with_ontology_columns(self):
        app = self.app()
        self.assertEqual(app.selectbox(key="project").value, "RPMM0003")
        self.assertFalse(app.checkbox(key="include_unlinked").value)
        self.assertEqual(len(app.tabs), 0)
        self.assertEqual(len(app.metric), 0)
        frame = self.column_frame(app)
        self.assertEqual(len(frame), 4)
        self.assertTrue(frame["Ontology ID"].astype(bool).all())
        self.assertEqual(set(frame.columns).intersection({"Column", "Ontology", "Ontology ID", "Source file"}),
                         {"Column", "Ontology", "Ontology ID", "Source file"})

    def test_concentration_counts_sum_context_mappings(self):
        app = self.app()
        values = self.inspect(app, "RPMM0068", "conc")
        self.assertEqual(self.counts(values), {"0": 24, "200": 24})
        self.assertEqual(values.columns.tolist(), ["Value", "Occurrences"])
        app.text_input(key="value_search").set_value("200").run()
        self.healthy(app)
        self.assertEqual(self.counts(self.value_frame(app)), {"200": 24})

    def test_same_header_in_two_source_files_remains_two_columns(self):
        app = self.app()
        acetate = self.inspect(app, "RPMM0027", "treatment", "AcetateTreatment")
        self.assertEqual(self.counts(acetate), {"AcOH": 5, "Water": 5})
        frame = self.column_frame(app)
        treatment_rows = frame[frame["Column"] == "treatment"]
        self.assertEqual(len(treatment_rows), 2)
        self.assertEqual(treatment_rows["Source file"].nunique(), 2)
        drought = self.inspect(app, "RPMM0027", "treatment", "DroughtTreatment")
        self.assertEqual(self.counts(drought), {"dry": 30, "water": 39})

    def test_null_genotype_is_a_recorded_literal_with_42_occurrences(self):
        app = self.app()
        values = self.inspect(app, "RPMM0063", "Genotype2")
        self.assertEqual(self.counts(values)["null"], 42)
        self.assertEqual(values.columns.tolist(), ["Value", "Occurrences"])

    def test_changing_project_resets_a_foreign_selected_column(self):
        app = self.app()
        self.inspect(app, "RPMM0068", "conc")
        app.selectbox(key="project").set_value("RPMM0027").run()
        self.healthy(app)
        self.assertTrue(app.selectbox(key="column").value.startswith("RPMM0027|"))
        expected = filter_columns(self.data.columns, project_ids=["RPMM0027"], mapped_only=True)
        self.assertEqual(len(self.column_frame(app)), len(expected))

    def test_zero_concept_project_can_reveal_its_unlinked_columns(self):
        app = self.app()
        app.selectbox(key="project").set_value("RPMM0012").run()
        self.healthy(app)
        self.assertEqual(len(self.value_frame(app)), 0)
        self.assertTrue(any("no ontology" in message.value.lower() for message in app.info))
        app.checkbox(key="include_unlinked").set_value(True).run()
        self.healthy(app)
        expected = filter_columns(self.data.columns, project_ids=["RPMM0012"])
        frame = self.column_frame(app)
        self.assertEqual(len(frame), len(expected))
        self.assertEqual(frame["Source file"].nunique(), 5)
        self.assertTrue(app.selectbox(key="column").value.startswith("RPMM0012|"))

    def test_hiding_unlinked_fields_clears_the_previous_value_search(self):
        app = self.app()
        app.checkbox(key="include_unlinked").set_value(True).run()
        app.selectbox(key="column").set_value(self.key("RPMM0003", "gcid")).run()
        app.text_input(key="value_search").set_value("FD_01").run()
        app.checkbox(key="include_unlinked").set_value(False).run()
        self.healthy(app)
        self.assertEqual(app.text_input(key="value_search").value, "")
        self.assertNotEqual(app.selectbox(key="column").value, self.key("RPMM0003", "gcid"))
        self.assertGreater(len(self.value_frame(app)), 0)

    def test_empty_value_search_does_not_show_unrelated_values(self):
        app = self.app()
        self.inspect(app, "RPMM0068", "conc")
        key = app.selectbox(key="column").value
        app.text_input(key="value_search").set_value("no_such_value_194857").run()
        self.healthy(app)
        self.assertEqual(len(self.value_frame(app)), 0)
        self.assertEqual(app.selectbox(key="column").value, key)
        self.assertTrue(any("no values" in message.value.lower() for message in app.info))

    def test_local_genotype_values_remain_available_under_a_schema_concept(self):
        app = self.app()
        values = self.inspect(app, "RPMM0068", "genotype")
        self.assertEqual(int(values["Occurrences"].sum()), 48)
        self.assertTrue({"Col-0", "ahk2ahk3"}.intersection(set(values["Value"])))
        self.assertEqual(values.columns.tolist(), ["Value", "Occurrences"])


if __name__ == "__main__":
    unittest.main()
