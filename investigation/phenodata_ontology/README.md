# RPMM Phenodata ontology annotations

This annotator processes files named exactly `phenodata.csv`, ignoring case,
under `projects/RPMM####/data`. Other `pheno*.csv` variants are outside
this run. Original CSVs are read only.

Run from any directory with Python 3.9 or later; no third-party packages or
network access are needed for an ordinary annotation run:

```sh
python3 investigation/annotate_phenodata.py
PYTHONPATH=investigation python3 -m phenodata_ontology.validate_export
python3 -m unittest investigation.phenodata_ontology.test_annotation
```

The commands above are shown relative to the repository root. Absolute script
paths also work. If the Python installation redirects bytecode caches outside
the workspace, set `PYTHONDONTWRITEBYTECODE=1`.

For a pilot, specify a separate output directory. A project filter without an
explicit output directory is rejected so it cannot overwrite the full export:

```sh
python3 investigation/annotate_phenodata.py RPMM0068 RPMM0027 RPMM0003 RPMM0017 --output /private/tmp/rpmm-phenodata-pilot
PYTHONPATH=investigation python3 -m phenodata_ontology.validate_export --output /private/tmp/rpmm-phenodata-pilot
```

## Output tables

The aggregate output is `reports/phenodata/`. Each project's filtered tables
are also generated in `projects/RPMM####/phenodata/`. `provenance.json` in that
folder links to the aggregate QA report and records shared rule/catalog hashes.
The validator checks that all 31 project exports exactly match their aggregate
subsets. Global external-source issues belong only to the aggregate curation table.

| File | Purpose |
| --- | --- |
| `inventory.tsv` | Source paths, SHA256, published and actual MD5, download URL, file IRI, row counts |
| `column_mappings.tsv` | One record per source column, its attribute category, type, unit and semantic rule |
| `value_mappings.tsv` | Context-scoped unique values and their annotations, expansion count and source records |
| `cell_annotations.tsv` | One record per source cell, including empty cells, typed values and ontology annotations |
| `row_annotations.tsv` | Source record identity, standard/reference classification and local sample/assay identifiers |
| `curation_queue.tsv` | Column definitions, unresolved values, contradictions and external-source problems |
| `coverage_by_project.tsv` | Cell status and ID assignment counts for each project |
| `coverage_by_role.tsv` | Value-level and cell-level coverage for each semantic role |
| `qa_report.json` | Census, status counts, source and rule hashes, runtime checks |
| `export_validation.json` | Independent read-back validation of the saved TSVs |

TSVs use standard CSV quoting with a tab delimiter. Read them with a CSV/TSV
parser: some original headers contain embedded newlines. Structured `assertions`,
`candidate_terms` and `source_records` fields contain JSON.

`source_record` is one-based and includes the header as record 1. It is distinct
from the physical `source_line_start` / `source_line_end`; a quoted header or
value can span several lines. `column_index` is one-based. The source file hash,
record number and column index identify a cell even when a local assay ID is
missing or reused. Local IDs do not assert biological-replicate identity.

`source_path` remains the stable logical key `RPMM####/...`; it is not a
repository-relative physical path. `inventory.tsv` now also has
`repository_path`, e.g. `projects/RPMM0003/data/AnalysedDataset/phenodata.csv`.
Logical keys and annotation evidence are unchanged by relocation, preserving
all original cell, row and mapping IDs. Custom `--output` locations write
project exports inside `OUTPUT/projects/` unless `--project-output-root` is set.

## Annotation semantics

Column terms describe attribute categories, not RDF predicates. A column's term
and a cell's term are different assertions. Numeric values, units, booleans,
dates, identifiers and text have explicit types; an ontology class is not
invented for every literal.

| Status | Meaning |
| --- | --- |
| `accepted` | A supported typed value or verified ontology annotation; this is rule-based acceptance, not a claim of external human review |
| `local` | A retained identifier, group code, strain name or descriptive literal |
| `needs_review` | Meaning, source context, unit or mapping remains unresolved |
| `missing` | Empty or explicitly unavailable data; original token is retained |
| `not_applicable` | Explicit standard-mixture marker in an organism/anatomy field |

Assigned ontology concepts use `term_id`; quantitative units use `unit_id`.
The normalized numeric value is separate from its unit and temporal reference
event. Signed flowering offsets are timepoints rather than negative durations.
`candidate_term_id`, `candidate_unit_id` and `candidate_terms` are suggestions,
not assignments. Component `assertions` also carry their own status. Review
candidates in cells never retain an assigned `term_id` or `unit_id`. Column
schemas may retain a verified attribute category or explicitly stated unit
while the column's precise meaning or measurement basis still needs review.

Values are scoped to file, column and context, including the row's standard or
reference classification. A project species is never copied onto all rows.
`Arabidopsis` is not narrowed to *A. thaliana* without explicit title or same-row
Col-0 evidence. Allium epithets use the same-row genus. Genotype `null` is retained
for review instead of automatically becoming missing. `Col-0` is normalized as a
local line name; stock accessions are not inferred. Drug annotations explicitly
describe the chemical agent, separately from the exposure condition.

## Vocabulary and evidence

`data/ontology_terms.json` is the checked-in vocabulary snapshot.
`data/ontology_sources.json` records official source URLs, release versions,
immutable GitHub commits, hashes and the original targeted OLS API responses.
The catalog checks its own hashes when loaded. All assigned IDs must be current
and verified in this snapshot. Obsolete `OBI:0100026` remains in the catalog as
historical evidence; new organism-category annotations use its verified
replacement `COB:0000022`.

Updates are explicit and require network access:

```sh
PYTHONPATH=investigation python3 -m phenodata_ontology.catalog --refresh
python3 investigation/annotate_phenodata.py --refresh-rdf
```

Changing a vocabulary or a rule requires regenerating and validating the export.
The export records the hashes of the vocabulary, mapping code, project-title
metadata, download manifest and RDF discovery evidence.

PMM file IRIs and download URLs come from the existing manifest. The first run
could not retrieve live PMM SPARQL results (query and connection checks timed
out); `data/rdf_discovery.json` records the query, attempt and failure. No
sample-level RDF annotation has been reused. `--refresh-rdf` retries file-parent
discovery only: a successful discovery still requires verified experiment,
measurement and sample joins before source RDF annotations can be inherited.

Two inputs, RPMM0037 and RPMM0068 LC, already had a published-MD5 mismatch in the
original download manifest. Their actual SHA256 identifies the source used for
annotation; the discrepancy is retained in inventory and curation outputs.

## Review workflow

Review `curation_queue.tsv` by `issue_code`, then follow its file/column/value
coordinates or `mapping_id` into the corresponding table. The RPMM0017 mixture
ratios conflict with their notes in 30 records; neither interpretation is chosen
automatically. Other examples include undocumented abbreviations, composite
descriptions, unspecified units and sampling-time reference events.

Confirm a rule from the source protocol/publication before updating
`columns.py` or `values.py`. Retain source text, cite the evidence in the rule,
add a targeted regression test, regenerate the annotations, then run the
read-back validator. Coverage is reported for both unique values and expanded
cells, separately from retained local literals. A high local-literal count is
expected because assay IDs and filenames are not ontology classes.
