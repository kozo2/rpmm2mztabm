# RPMM ontology columns — Streamlit

A focused Python browser for checking each RPMM project's ontology columns and
the values they contain.

1. Choose a project.
2. Check its ontology-column list and select a column.
3. Read the original values and occurrence counts, or download them as CSV.

Enable **Include columns without a concept ID** to inspect local metadata and
projects whose columns have not been linked to ontology concepts. Source files
remain separate when they contain columns with the same name. Column details
are available in a collapsed section.

## Run

Requires Python 3.12 or later. From `dashboard/`:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/streamlit run app.py
```

Open the Local URL printed by Streamlit (normally http://localhost:8510).
No Node.js or TypeScript build is required.

## Data

Inside the rpmm2mztabm repository the app reads the authoritative reports:

- `reports/phenodata/column_mappings.tsv`: source columns and schema concepts
- `reports/phenodata/value_mappings.tsv`: raw and normalized values, annotations and cell counts
- `reports/phenodata/inventory.tsv`: source file provenance
- `investigation/projects.json`: project titles

A bundled `snapshot/` supports a standalone dashboard checkout. To refresh it
from the scientific reports:

```sh
python3 scripts/build-data.py
python3 scripts/test-data.py
```

The generator preserves original source keys and numeric strings, records source
SHA256 values, and produces identical bytes for unchanged reports. See
[scripts/DATA.md](scripts/DATA.md) for the data contract.

Column identity includes project, source file and column index. Column concepts
are schema categories. The values table groups exact original strings and sums
their source-cell counts, including repeated mappings with different row context.
The underlying reports and snapshot retain all mapping context and candidates.

The complete dataset has 31 projects, 39 phenodata files, 391 column instances
and 11,784 contextual value mappings covering 26,997 cells. Three projects with
no column concept IDs still have inspectable metadata.
