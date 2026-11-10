# Plan: Convert `projects/<RPMM_ID>/data` (RPMM) to mzTab-M 2.1 with RmzTabM

Status: **plan, feasibility verified** — the reference implementation
[`mztabm/proto_rpmm0003.R`](mztabm/proto_rpmm0003.R) already produces a file that passes
the official jmztabm validator with **0 errors / 0 warnings / 0 messages** and round-trips
with identical abundance values.


Current repository layout (relocated 2026-11-10): source files are in
`projects/RPMM####/data/`; project artifacts are adjacent in `tables/`,
`phenodata/`, and `mztabm/`. Cross-project reports live in `reports/`. The
prototype and conversion source code remain in `investigation/`. Source keys
such as `RPMM0003/AnalysedDataset/MetMat_abs.txt` remain stable logical identifiers;
they resolve physically to `projects/RPMM0003/data/AnalysedDataset/MetMat_abs.txt`.

Re-run the existing prototype from any working directory with:

```sh
Rscript /absolute/path/to/rpmm2mztabm/investigation/mztabm/proto_rpmm0003.R --output /private/tmp/RPMM0003_MetMat_abs.mztab
```

It accepts `--projects-root` and `--validator` (or `RPMM_MZTABM_VALIDATOR`). The
reported R/package feasibility below records the original investigation; install
those dependencies before re-running it. No mzTab-M output is regenerated merely
by relocating the files.

---

## 1. Goal and scope

Convert the downloaded RPMM (RIKEN Plant Metabolome Metadatabase) analysed datasets in
`projects/<RPMM_ID>/data/` into valid **mzTab-M 2.1.0** files using the
**RmzTabM** package in the local R environment.

In scope:

- Every per-project **quantification matrix** (peak-height / normalised data matrix),
  one mzTab-M file per matrix.
- Linked **sample metadata** (`phenodata*.csv`) → MTD samples/assays/study variables.
- Linked **NIST MSP annotation files** (`*_NIST.txt`) → SME identification evidence +
  feature m/z source, plus opt_ columns (RI, annotation grade, chemical name).
- `ALKANE_SERIE.txt` (RI calibration) → recorded as an opt_/comment provenance note.

Out of scope (kept as companion files, not mzTab-M content):

- `*limma*_Diff_res*` differential statistics (mzTab-M has no stats section).
- `.RData`, `.DS_Store`, publication supplements not used as quant matrices.

Open decision (default chosen): convert **all** quant matrices per project (both
`MetMat_abs`-style raw matrices and normalised ones), so the normalisation state is
explicit in the filename, e.g. `RPMM0003_MetMat_abs.mztab` and `RPMM0003_CRMN_normalized.mztab`.

---

## 2. Verified environment facts

| Item | Value |
|---|---|
| R | 4.6.1 (`/usr/local/bin/R`) |
| RmzTabM | 0.99.3 (installed; loaded from `RmzTabM`) |
| SummarizedExperiment | available |
| data.table / jsonlite | available |
| readxl / openxlsx | **NOT installed** — needed for `.xlsx`/`.xls` matrices (RPMM0032 is xlsx-only) |
| Official validator | jmztabm-cli 2.1.0-SNAPSHOT binary already cached at `~/Library/Caches/org.R-project.R/R/BiocFileCache/8be021b8c274_jmztabm-validator-aarch64-macos`; `writeMzTabM(..., validate = TRUE)` and `readMzTabM()` call it automatically (needs Java — present, it ran during prototyping) |

---

## 3. Data landscape (census)

- **47 project folders** `RPMM0001`–`RPMM0068`, each with one or more assay subfolders
  (mostly `AnalysedDataset/`; some projects split by technique:
  `RPMM0004/{ITEX,ITEX_and_SPME}`, `RPMM0018/Year2/{gc012,gc2}`, `RPMM0061/{long,short}`,
  `RPMM0068/{GC,LC,Lipid}`, `RPMM0026/{LipophilicProfile,PolarMetabolites}`, …).
- ~282 files total. Roles (by filename census):
  - **Quant matrices** (~60): `MetMat_abs.txt` (25), `CRMN_normalized.txt` (19+4 variants),
    `one_normalized.txt` (8+2), plus ~15 one-off shapes (LC peak lists, Kusano
    named-metabolite matrices, xlsx/xls matrices).
  - **Sample metadata** (~60): `phenodata*.csv`, `Phenodata*.csv`, `pheno.csv`, …
    (column names vary; sample-ID column usually `gcid`, matching matrix column names).
  - **NIST MSP annotations** (~30): `*NIST.txt` — MSP blocks with `Synon:` fields
    (`RI`, `RT(s)`, `Annotated`, `Protocol`, `Species`, `PI`, `Date`, `MCR-Reference`)
    and EI peak lists; `MCR-Reference` matches the matrix `Win` column tokens.
  - **Alkane series** (~51): `ALKANE_SERIE.txt`.
  - **limma results** (~10), **Excel supplements** (~14).
- Technique mix: predominantly GC-MS (EI, retention time + RI, no precursor m/z);
  LC-MS and lipidomics variants exist (e.g. `RPMM0068/AnalysedDataset/LC/94_CYTresult_TM_RN_pos.txt`
  has `Mass`, `Ret(umin)`, `Annotation`, `Compound` columns).

---

## 4. mzTab-M mapping design (per project → one file per matrix)

Profile: **M+S+F+E** (the only profile in RmzTabM 0.99.3 that allows SML+SMF together).
Sections:

| mzTab-M section | Source | Notes |
|---|---|---|
| MTD `mzTab-ID` | RPMM ID + matrix name | e.g. `RPMM0003_MetMat_abs` |
| MTD samples | `phenodata` rows | species as CV term (e.g. `[NCBITaxon, NCBITaxon:4081, Solanum lycopersicum, ]`); per-project species curated in a config table |
| MTD ms_run[1..n] | per sample | location = source-file URL + `#<sample_id>` fragment (unique basename, honest target); `scan_polarity` plain `"positive"/"negative"` |
| MTD assays | 1:1 with samples | `assay_cols` → `abundance_assay[n]` |
| MTD study variables | phenodata group columns (e.g. `class1`, `fruitvsleaf`) | via `groups=` of `MzTabM(se)`; drives `abundance_study_variable[n]` in SML |
| MTD database | NIST | `setMtdDatabase(name="[,,NIST,]", prefix="NIST", ...)` |
| MTD `id_confidence_measure[1]` | `[MS, MS:1002896, compound identification confidence level, ]` | added by direct MTD row edit (see §6) |
| SME | one evidence row per feature | `chemical_name` from matrix annotation column; `database_identifier = "NIST"` for annotated rows, `"null"` otherwise; `identification_method = "[,,NIST MS library search,]"`; `spectra_ref = "ms_run[1]:index=0"`; `rank=1` |
| SMF | one feature row per matrix row | `retention_time_in_seconds` from RT column; `SME_ID_REFS` → SME_ID integers; abundances → `abundance_assay[n]`; opt_ columns: RI, annotation grade, chemical name, window token |
| SML | 1:1 with SMF (features are already collapsed) | `chemical_name`, `reliability`, `database_identifier`; + `abundance_study_variable[n]` via `smlAddStudyVariableColumns()` |

Key semantic decisions (documented in-file via opt_ columns):

1. **Feature m/z (SMF/SME `exp_mass_to_charge`, mandatory + non-nullable):** GC-MS rows have
   no precursor m/z, so use the **base-peak m/z of the deconvoluted EI spectrum** from the
   project's NIST MSP file, matched via the `Win`/`MCR-Reference` token (RPMM0003 coverage:
   34/39 directly; remainder via nearest-RI fallback). Documented per row with
   `opt_global_mz_source = "base peak m/z of deconvoluted EI spectrum"`.
   For LC-MS matrices with a real `Mass` column, use that column instead.
2. **SME `charge` = 1** (EI product ions / singly-charged features); SMF `charge` stays `null` (nullable there).
3. **SME `theoretical_mass_to_charge`** = library base-peak m/z (mandatory column; value is the library-side expectation, flagged by `opt_global_mz_source`).
4. **Reliability mapping** `IIS_I_A_or_U` → mzTab-M levels: `A2 → 2`, `I → 3`, `U → 4`;
   original grade preserved in `opt_global_annotation_grade`.
5. Empty/`NA` abundances are written as `null` (allowed in abundance columns).

---

## 5. RmzTabM 0.99.3 quirks catalogued (all hit & solved during prototyping)

These must be handled by the converter lib — they are why a naive script fails:

1. `MzTabM(se)` (SummarizedExperiment path) emits **only MTD** unless a non-empty
   `smfCols.` is passed; the SE path also can't produce SML/SME — assemble manually.
2. Supported profiles are only `M, M+S, M+F, M+F+E, M+S+F+E` — SML+SMF requires **M+S+F+E**
   (i.e. an SME section must be supplied too).
3. `smeCreate()` emits an `adduct_ions` column in SME; the spec/validator requires
   **`adduct_ion`** (singular) → rename post-hoc.
4. SMF and SME `exp_mass_to_charge` (and SME `charge`, `theoretical_mass_to_charge`) are
   **mandatory in the header AND non-nullable** in cells → must carry real values (§4.1–3).
5. SME **list columns must have consistent item counts** (`chemical_name` vs
   `database_identifier`) → `database_identifier = "NIST"` exactly on annotated rows.
6. `setMtdField()` **cannot create `id_confidence_measure[n]`** (it would emit
   `id_confidence[n]`) → rbind the MTD row directly + `mtdSort()`.
7. `smeCreate()` **requires** `exp_mass_to_charge`/`charge`/`theoretical_mass_to_charge`
   args even if you drop the columns afterwards.
8. Unique **ms_run locations** are needed for R-side round-trip row names → `#<sample>` fragments.
9. `scan_polarity` in colData must be plain `"positive"/"negative"` (not a CV string);
   `species`/`tissue` are CV strings.
10. Read-side: `readMzTabM()` leaves SMF rownames empty and
    `makeSummarizedExperimentFromMzTabM()` mangles `#` colnames via `make.names()`
    → QA must compare parsed sections directly (abundance equality), not via the SE helper.
11. Matrix header off-by-one risk: RPMM0003-style matrices have **7 annotation columns**
    before sample columns; do not assume a fixed count — detect by intersecting with phenodata IDs.

---

## 6. Pipeline architecture

```
investigation/mztabm/
  proto_rpmm0003.R        # DONE - reference implementation, all quirks handled
  lib/
    inventory.R           # classify every file: role, matrix shape, phenodata variant, MSP presence
    read_matrix.R          # readers per matrix shape (see §7 families)
    read_phenodata.R      # phenodata variants -> standard colData frame
    parse_msp.R           # NIST MSP parser: base peak, RI, RT, Synon fields, per-window table
    build_mztabm.R       # MTD/SMF/SME/SML assembly incl. all §5 workarounds
    qa.R                  # validator invocation + round-trip checks + QA-report rows
  projects.tsv             # per-project config: species CV, polarity, group columns,
                          #   matrix->phenodata pairing, overrides
  run_conversion.R         # batch driver; skips completed, collects failures
Output:
  projects/<RPMM_ID>/mztabm/<RPMM_ID>_<MatrixName>.mztab
  reports/mztabm/QA-report.tsv          # per-file validator status + counts
  reports/mztabm/conversion-errors.log  # per-file exception + reason
```

Per file pipeline: detect shape → read matrix + phenodata (+ MSP) → build SE for MTD →
set MTD database/id_confidence_measure → SME → SMF → SML (+ study variable columns) →
`writeMzTabM(validate = TRUE)` → QA re-read + abundance equality → log.

Failure policy: any file that cannot meet a mandatory field honestly is **skipped and
logged** (never fabricated silently); it lands on a curation list reviewed per §7 Phase 4.

---

## 7. Execution phases

**Phase 0 — done.** Feasibility prototype on RPMM0003 (`proto_rpmm0003.R`): validator-clean
M+S+F+E file; all API quirks documented above.

**Phase 1 — inventory & config (half day).**
- `inventory.R`: walk all 47 projects, classify every file, emit a conversion manifest
  (matrix file → phenodata file → MSP file → assay subfolder) with per-shape family tags.
- Curate `projects.tsv`: species CV term per project (from MSP `Synon:Species` +
  Metabolonote metadata), polarity, study-variable group columns, matrix/phenodata pairing
  for non-standard names.

**Phase 2 — GC-MS "MetMat family" converter (~1 day).** Covers the dominant shape
(`MetMat_abs`, `CRMN_normalized`, `one_normalized`, …: annotation columns + RI + Win +
sample columns, with phenodata + MSP siblings). Generalise `proto_rpmm0003.R` into
`lib/` modules and `run_conversion.R`. Expected coverage: ~30 projects / ~55 matrices.

**Phase 3 — LC-MS / lipidomics shapes (~1 day).** e.g. `RPMM0068/{LC,Lipid}`
(`Mass`, `Ret(umin)`, `Annotation`, `Compound` columns → real `exp_mass_to_charge`,
RT conversion µmin→s), `RPMM0032` xlsx (needs `readxl` installed first:
`install.packages("readxl")`), `RPMM0013`, `RPMM0026`.

**Phase 4 — stragglers & curation (~1 day).** Kusano named-metabolite matrices
(`RPMM0006`: no RT/RI columns → RT null, m/z from MSP where possible), `.xls` legacy files,
projects without MSP (fallback m/z policy review), anything Phase 2–3 logged as skipped.

**QA gate (every phase):** official validator 0 errors (warnings reviewed), re-read
abundance equality vs source matrix, sample count and feature count checks, QA-report
aggregates printed. Spot-check one file per phase manually in a text editor.

Estimated total: ~3–4 working days after plan approval.

---

## 8. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Validator evolves / strictness differences | pin to the cached jmztabm-cli binary; QA-report keeps per-file status |
| Fabricated-looking m/z for GC-MS features | §4.1 policy: base-peak m/z only, always flagged via `opt_global_mz_source`; alternative (omit) is impossible under 2.1 mandatory-column rules |
| Format drift across 47 projects | shape families + per-project overrides in `projects.tsv`; failures logged, never guessed |
| RmzTabM 0.99.3 bugs (§5) | workarounds centralised in `build_mztabm.R`; revisit if package updates (>= 1.0) fix them |
| xlsx/xls reading | install `readxl` (Phase 3); xls fallback via `libxls`/conversion script if needed |
| Large matrices (e.g. 175+ samples) | file-based writing; validator handles (RPMM0003 full = 183 cols OK structurally; only runtime scales) |

## 9. Definition of done

- Every convertible matrix (target: all ~60) has an mzTab-M 2.1.0 file under
  `projects/<RPMM_ID>/mztabm/`.
- `QA-report.tsv`: 0 validator errors per emitted file; abundance round-trip equality.
- Unconverted files explicitly listed with reasons in `conversion-errors.log`.
- One representative file per shape family manually reviewed.
