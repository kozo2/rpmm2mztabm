"""A focused project → ontology column → recorded values Streamlit browser."""

from pathlib import Path
import sys

import streamlit as st

DASHBOARD = Path(__file__).resolve().parent
sys.path.insert(0, str(DASHBOARD))

from atlas_data import get_values, load_data, to_csv

st.set_page_config(page_title="RPMM ontology columns", page_icon="🌱", layout="wide")


def source_fingerprint():
    reports = DASHBOARD.parent / "reports" / "phenodata"
    paths = [reports / name for name in ("column_mappings.tsv", "value_mappings.tsv", "inventory.tsv")]
    if all(path.is_file() for path in paths):
        paths.append(DASHBOARD.parent / "investigation" / "projects.json")
    else:
        paths = sorted((DASHBOARD / "snapshot").rglob("*.json"))
    return tuple((str(path), path.stat().st_mtime_ns, path.stat().st_size)
                 for path in paths if path.exists())


@st.cache_data(show_spinner="Loading columns and values…")
def cached_data(fingerprint):
    return load_data()


def reset_column():
    st.session_state.pop("column", None)
    st.session_state.pop("value_search", None)


st.title("RPMM ontology columns")
st.caption("Choose a project, check its ontology columns, and see the values they contain.")

try:
    data = cached_data(source_fingerprint())
except Exception as error:
    st.error(f"The column mappings could not be loaded: {error}")
    st.stop()

project_lookup = data.projects.set_index("id").to_dict("index")
project_ids = sorted(project_lookup)
project = st.selectbox("Project", project_ids, key="project", on_change=reset_column,
                       format_func=lambda value: f"{value} · {project_lookup[value]['title']}")

include_unlinked = st.checkbox("Include columns without a concept ID", key="include_unlinked",
                               on_change=reset_column,
                               help="Also show local identifiers, free text and unit-only columns.")
columns = data.columns.loc[data.columns["project"].eq(project)].copy()
if not include_unlinked:
    columns = columns.loc[columns["columnTermId"].ne("")]
columns = columns.sort_values(["sourcePath", "columnIndex"], kind="stable").reset_index(drop=True)

st.subheader("Columns")
if columns.empty:
    st.info("No ontology-linked columns are recorded for this project. Enable “Include columns without a concept ID” to inspect its other phenodata columns.")
    st.stop()

column_table = columns[["sourceColumn", "columnTermLabel", "columnTermId", "sourcePath"]].rename(
    columns={"sourceColumn": "Column", "columnTermLabel": "Ontology", "columnTermId": "Ontology ID", "sourcePath": "Source file"})
column_table["Source file"] = column_table["Source file"].str.removeprefix(project + "/")
st.dataframe(column_table, hide_index=True, width="stretch", height=min(330, 38 + len(column_table) * 35),
             column_config={"Column": st.column_config.TextColumn(width="small"),
                            "Ontology ID": st.column_config.TextColumn(width="small"),
                            "Source file": st.column_config.TextColumn(width="large")})
st.caption(f"{len(columns)} source columns across {columns['sourcePath'].nunique()} file(s). Ontology IDs describe column categories.")

st.subheader("Values")
column_lookup = columns.set_index("key").to_dict("index")
keys = columns["key"].tolist()
if st.session_state.get("column") not in keys:
    st.session_state["column"] = keys[0]


def column_label(key):
    item = column_lookup[key]
    return f"{item['sourceColumn']} · {item['sourcePath'].removeprefix(project + '/')}"


column_key = st.selectbox("Column", keys, key="column", format_func=column_label,
                          on_change=lambda: st.session_state.pop("value_search", None))
column = column_lookup[column_key]
values = get_values(data, column_keys=[column_key])
frequencies = values.groupby("rawValue", dropna=False, sort=False)["cellCount"].sum().reset_index()
frequencies = frequencies.sort_values(["cellCount", "rawValue"], ascending=[False, True], kind="stable")
frequencies = frequencies.rename(columns={"rawValue": "Value", "cellCount": "Occurrences"}).reset_index(drop=True)
search = st.text_input("Search values", key="value_search", placeholder="Find a recorded value…")
if search.strip():
    frequencies = frequencies.loc[frequencies["Value"].astype(str).str.contains(search.strip(), case=False, regex=False)]

if frequencies.empty:
    st.info("No values match this search.")
else:
    st.dataframe(frequencies, hide_index=True, width="stretch",
                 height=min(460, 38 + len(frequencies) * 35),
                 column_config={"Value": st.column_config.TextColumn(width="large"),
                                "Occurrences": st.column_config.NumberColumn(format="%d", width="small")})
    st.caption(f"{len(frequencies)} distinct values · {int(frequencies['Occurrences'].sum())} occurrences. Values are shown as recorded; blank entries are empty cells.")
    st.download_button("Download values CSV", to_csv(frequencies),
                       file_name=f"{project}-column-{column['columnIndex']}-values.csv",
                       mime="text/csv", key="download_values", icon=":material/download:")

with st.expander("Column details"):
    st.write(f"**Source:** {column['sourcePath']} · column {column['columnIndex']}")
    if column.get("columnTermId"):
        st.write(f"**Ontology:** {column['columnTermLabel']} · {column['columnTermId']}")
    if column.get("unitId"):
        st.write(f"**Declared unit:** {column.get('unitLabel', '')} · {column['unitId']}")
    st.write(f"**Mapping state:** {str(column['status']).replace('_', ' ')}")
    if column.get("reason"):
        st.caption(column["reason"])
