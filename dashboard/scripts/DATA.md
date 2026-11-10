# Dashboard data snapshot

The dashboard reads a deterministic snapshot generated from the repository's
annotation reports. It does not classify columns, infer ontology matches or
query live ontology services.

From the repository root:

```sh
python3 dashboard/scripts/build-data.py
python3 dashboard/scripts/test-data.py
```

Python 3.9 or later and the standard library are sufficient. Explicit
`--reports-root`, `--projects-file` and `--output-dir` options are available.
Generation does not edit the original reports or annotations.

## Input and output

- `reports/phenodata/column_mappings.tsv` defines the available source columns.
- `reports/phenodata/value_mappings.tsv` supplies recorded values, context,
  annotations and weighted cell counts.
- `reports/phenodata/inventory.tsv` supplies file provenance and record counts.
- `investigation/projects.json` supplies project titles, with HTML removed.

`snapshot/explorer.json` contains the complete column index and project,
source-file and ontology-group metadata. Value mappings are split into
`snapshot/values/RPMM####.json` so the app can use a committed fallback when the
original repository reports are not available.
`snapshot/column_mappings.tsv` is an exact copy for download. Each input file's
SHA256 is included in the index, and the combined source fingerprint has no
timestamp. Regenerating unchanged sources produces identical output bytes.

The Streamlit app's `atlas_data.py` prefers the live repository reports. If those
files are absent, it reads this committed snapshot instead. A partly missing or
invalid report set is reported as an error rather than silently hidden behind a
snapshot. Pure filtering, matrix, aggregation and CSV helpers do not import
Streamlit and can be tested independently.

The checked-in snapshot contains 31 projects, 39 CSV files, 391 source column
instances, 165 original header spellings, 112 standardized field names, 11,784
context-scoped value mappings and 26,997 cells. There are 143 source columns with
an ontology category and 14 distinct category IDs. Of those 143 columns, 100
have status `accepted`, 23 `local` and 20 `needs_review`; category presence alone
does not mean that every detail of a column's meaning is resolved.

## Identity and counts

The column key is `project|source_path|column_index`. Source paths and original
headers are preserved, including embedded newlines and trailing whitespace.
Neither the header nor an ontology term alone uniquely identifies a column.

`valueCount` counts context-scoped mapping records. `distinctRawValues` counts
distinct raw strings, including the empty string. `cellCount` sums the original
`cell_count`, so repeated strings in different sample contexts are counted
correctly. `valueStatusCounts` counts mappings; `cellStatusCounts` is weighted by
cell counts. `assignedTermCount` and `assignedUnitCount` count distinct cell-level
IDs, while `termCellCount` and `unitCellCount` count affected cells.

The canonical group key uses `term:<column ontology ID>` when the source column
has a category and otherwise `field:<standard_name>`. Groups retain all of their
standardized field names: one generic category such as mass can describe several
different recorded quantities. The UI can compare either grouping without
rewriting source identities.

## Scientific interpretation

Column categories and cell-value terms describe different things. Candidate
terms/units, structured assertions, row kinds, time reference events and review
reasons are preserved separately from assigned IDs. Numeric annotations remain
decimal strings, avoiding a float conversion. Genotype `null` is retained and
is not turned into a missing value.

RPMM0012, RPMM0023 and RPMM0032 have valid phenodata files and local columns but
no column ontology category IDs in this snapshot. A view filtered to mapped
categories should not imply that these projects lack metadata.
