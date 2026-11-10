#!/usr/bin/env python3
"""export_tables.py - per-matrix export of RPMM quantification data.

For every quantification matrix under projects/<RPMM_ID>/data, write:

  projects/<RPMM_ID>/tables/<RPMM_ID>_<MatrixName>.toml
      sample metadata for the mzTab-M MTD table: project title, source file,
      organism, instrument hints, one entry per sample (assay, ms_run
      location, phenodata study variables).

  projects/<RPMM_ID>/tables/<RPMM_ID>_<MatrixName>.parquet
      per-metabolite quantification table for the mzTab-M SML table:
      annotation columns (feature id, chemical name, grade, RI/RT, m/z when
      present) + one abundance_<sample> column per sample.

Fail-and-log policy: anything that cannot be determined from the data files
(unknown organism, unmatched phenodata ids, unclassifiable columns) is
recorded in reports/tables/_curation.tsv - never invented.

Usage: python3 investigation/export_tables.py [RPMM_ID ...]
"""

import argparse
import csv
import datetime
import math
import os
import re
import sys
from pathlib import Path
import tomli_w

import pyarrow as pa
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = REPO / "projects"
OUT = REPO / "projects"
REPORTS = REPO / "reports" / "tables"
PROJECTS_JSON = HERE / "projects.json"
GENERATED_ON = datetime.date.today().isoformat()


def project_data_dir(project):
    """Physical source directory; project identity is independent of its location."""
    return os.path.join(DATA, project, "data")


def logical_path(path):
    """Stable source key RPMM####/... (without the physical data/ container)."""
    parts = Path(path).resolve().relative_to(Path(DATA).resolve()).parts
    if len(parts) < 3 or parts[1] != "data":
        raise ValueError(f"not a project data path: {path}")
    return Path(parts[0], *parts[2:]).as_posix()


def physical_path(logical):
    """Resolve a recorded logical key under the selected projects root."""
    parts = Path(logical).parts
    if len(parts) < 2 or not re.fullmatch(r"RPMM\d+", parts[0]):
        raise ValueError(f"invalid logical project source key: {logical}")
    target = (Path(DATA) / parts[0] / "data" / Path(*parts[1:])).resolve()
    target.relative_to(Path(project_data_dir(parts[0])).resolve())
    return str(target)


def output_reference(path):
    """Portable repository-relative reference for ordinary output locations."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(REPO).as_posix()
    except ValueError:
        return str(resolved)

# ---------------------------------------------------------------------------
# Text helpers (mirrors investigation/mztabm/lib/read_matrix.R rules)
# ---------------------------------------------------------------------------

META_NAMES = {
    "numbering_for_matrix", "numbering_for_unique_peaks", "uniquepeak",
    "representative_peak_selected", "marked_peaks", "iis_i_a_or_u",
    "i_a_m_or_p", "iamp", "final_annotation", "ri", "ri_", "win", "tag",
    "lineno", "query", "known", "cas", "preferred", "action", "peak_nr",
    "peak", "scan_nr", "ret_umin", "ri_", "m_z", "mz", "mass", "mass_scan",
    "deconv", "annotation", "compound", "check", "mstag", "anplus", "parent",
    "child", "rt", "rt_min", "rt_s", "retention", "retention_time",
    "retention_time_min", "inchi", "smiles", "formula", "database", "name",
    "metabolites", "synonym", "biological_name", "num", "no", "index",
    "sample_id", "sample", "peak_name", "list", "sample_type", "kinds",
    "project_id", "authentic_sample", "cx", "x",
}
META_PREFIXES = ("fdata", "numbering", "marked", "iamp", "synon", "opt_")

HEADER_ANCHORS = {
    "win", "ri", "final annotation", "numbering for matrix", "iis_i_a_or_u",
    "ret(umin)", "peak nr", "annotation", "metabolites", "synonym", "name",
    "m_z", "m/z", "mass", "retention", "retention time", "biological name",
}


def norm_name(x):
    return re.sub(r"[^a-z0-9]+", "_", (x or "").strip().lower()).strip("_")


def is_meta_name(name):
    n = norm_name(name)
    return n in META_NAMES or n.startswith(META_PREFIXES)


def to_num(v):
    """Parse a cell to float, else NaN. Handles str/float cells uniformly."""
    if v is None:
        return math.nan
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v)
    s = str(v).strip().strip('"')
    if s in ("", "NA", "NaN", "null", "NULL", "#N/A", "nan"):
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def fmt_cell(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return str(int(v)) if v.is_integer() else repr(v)
    return str(v).strip()


def guess_delim(lines1):
    head = ""
    for ln in lines1:
        if ln.strip():
            head = ln
            break
    return "\t" if head.count("\t") >= head.count(",") else ","


def find_header_row(rows, max_scan=25):
    for i, row in enumerate(rows[:max_scan]):
        cells = [norm_name(c) for c in row]
        if sum(bool(c) for c in cells) < 3:
            continue
        low = [c.lower() for c in (x.strip() for x in row if x)]
        if any(c in HEADER_ANCHORS for c in low) or any(
                re.match(r"^ret\(|^numbering|^peak nr|^m/z", c) for c in low):
            return i
    for i, row in enumerate(rows[:max_scan]):
        if sum(bool(str(c).strip()) for c in row) >= 3:
            return i
    return 0


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_rows(path):
    """Return (rows: list[list[str]], format). Text = detected delimiter,
    Excel = first sheet with a locatable header."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xlsm"):
        import openpyxl
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
        ws = wb[wb.sheetnames[0]]
        rows = [[fmt_cell(c) for c in r] for r in ws.iter_rows(values_only=True)]
        wb.close()
        return trim_rows(rows), "xlsx"
    if ext == ".xls":
        import xlrd
        wb = xlrd.open_workbook(path)
        sh = wb.sheet_by_index(0)
        rows = [[fmt_cell(sh.cell_value(r, c)) for c in range(sh.ncols)]
                for r in range(sh.nrows)]
        return trim_rows(rows), "xls"
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    d = guess_delim(lines)
    rows = [[c.strip().strip('"') for c in ln.split(d)]
            for ln in lines if ln.strip()]
    return trim_rows(rows), "tsv" if d == "\t" else "csv"


def trim_rows(rows):
    width = max((len(r) for r in rows), default=0)
    out = []
    for r in rows:
        r = list(r) + [""] * (width - len(r))
        out.append([fmt_cell(c) for c in r])
    return out


# ---------------------------------------------------------------------------
# Matrix / phenodata classification
# ---------------------------------------------------------------------------

def read_matrix(path):
    rows, fmt = read_rows(path)
    hdr_i = find_header_row(rows)
    header = [c.strip() for c in rows[hdr_i]]
    body = rows[hdr_i + 1:]
    if not body:
        raise ValueError(f"no data rows: {path}")

    meta_idx, sample_idx, unknown_idx = [], [], []
    for j, name in enumerate(header):
        if not name:
            continue
        vals = [r[j] for r in body[:60]]
        numeric = sum(1 for v in vals if not math.isnan(to_num(v)))
        frac = numeric / max(1, len(vals))
        if is_meta_name(name):
            meta_idx.append(j)
        elif frac >= 0.5:
            sample_idx.append(j)
        else:
            unknown_idx.append(j)      # text, unrecognized -> not a sample
    if not sample_idx:
        raise ValueError(f"no sample columns recognized: {path}")

    meta = {}
    for j in meta_idx + unknown_idx:
        name = header[j] or f"col{j}"
        meta[norm_name(name)] = [r[j].strip() for r in body]
    samples = [header[j] for j in sample_idx]
    abund = {}
    for j, s in zip(sample_idx, samples):
        abund[s] = [to_num(r[j]) for r in body]
    n = len(body)
    return {
        "path": path, "format": fmt, "header": header, "n_rows": n,
        "meta": meta, "samples": samples, "abund": abund,
        "meta_columns": [header[j] for j in meta_idx],
        "unknown_columns": [header[j] for j in unknown_idx],
    }


def read_phenodata(path):
    rows, _ = read_rows(path)
    hdr_i = find_header_row(rows)
    header = [c.strip() for c in rows[hdr_i]]
    body = rows[hdr_i + 1:]
    cols = {header[j] or f"col{j}": [r[j].strip() for r in body]
            for j in range(len(header))}
    return cols


def align_phenodata(pcols, samples):
    """Pick the sample-id column by max overlap with matrix sample names."""
    best, best_hits = None, -1
    for name, vals in pcols.items():
        hits = sum(1 for s in samples if s in vals)
        bonus = 3 if norm_name(name) in (
            "gcid", "sample_id", "sample", "samplename", "sample_name",
            "id", "filename") else 0
        if hits + bonus > best_hits:
            best, best_hits = name, hits + bonus
    if best is None:
        return None, {}
    ids = pcols[best]
    aligned = {}
    for s in samples:
        if s in ids:
            i = ids.index(s)
            aligned[s] = {k: v[i] for k, v in pcols.items() if i < len(v)}
    return best, aligned


# ---------------------------------------------------------------------------
# Species resolution (curated map; unknown -> curation, never guessed)
# ---------------------------------------------------------------------------

SPECIES_MAP = {
    "arabidopsis": ("NCBITaxon:3702", "Arabidopsis thaliana"),
    "arabi": ("NCBITaxon:3702", "Arabidopsis thaliana"),
    "tomato": ("NCBITaxon:4081", "Solanum lycopersicum"),
    "rice": ("NCBITaxon:4530", "Oryza sativa"),
    "ine": ("NCBITaxon:4530", "Oryza sativa"),
    "lettuce": ("NCBITaxon:4236", "Lactuca sativa"),
    "sunnylettuce": ("NCBITaxon:4236", "Lactuca sativa"),
    "cassava": ("NCBITaxon:3983", "Manihot esculenta"),
    "chickpea": ("NCBITaxon:3827", "Cicer arietinum"),
    "soybean": ("NCBITaxon:3847", "Glycine max"),
    "fermented soybean": ("NCBITaxon:3847", "Glycine max"),
    "soy": ("NCBITaxon:3847", "Glycine max"),
    "natane": ("NCBITaxon:3708", "Brassica napus"),
    "brassica napus": ("NCBITaxon:3708", "Brassica napus"),
    "brassica rapa": ("NCBITaxon:3711", "Brassica rapa"),
}


def resolve_species(raw):
    v = re.sub(r"\s+", " ", (raw or "").strip().lower())
    v = re.sub(r"\s*\(.*$", "", v)
    v = re.sub(r"\s+(leaf|fruit|root|seed|stem|nodule).*$", "", v).strip()
    if not v:
        return None
    if v in SPECIES_MAP:
        return SPECIES_MAP[v]
    for k, cv in SPECIES_MAP.items():
        if v.startswith(k):
            return cv
    return None


def find_species_raw(matrix_dir, project_id, pdata):
    """Raw species string: MSP Synon:Species line, then phenodata columns."""
    for root, _dirs, files in os.walk(project_data_dir(project_id)):
        for fn in sorted(files):
            if "nist" in fn.lower() and fn.lower().endswith(".txt"):
                try:
                    with open(os.path.join(root, fn), encoding="utf-8",
                              errors="replace") as f:
                        for line in f:
                            if line.startswith("Synon:Species:"):
                                v = line.split(":", 2)[2].strip()
                                if v:
                                    return v, f"msp:{fn}"
                except OSError:
                    pass
    for col, vals in (pdata or {}).items():
        if re.search(r"species|binomial", col, re.I):
            vals = [v for v in vals if v]
            if vals:
                return vals[0], f"phenodata:{col}"
    return None, None


# ---------------------------------------------------------------------------
# Inventory
# ---------------------------------------------------------------------------

SKIP_RE = re.compile(r"\.DS_Store|\.rdata$|\.rda$", re.I)
COMPANION_RE = re.compile(r"limma|diff_res|_res\.", re.I)


def inventory():
    """project -> list of matrix files (de-duplicated txt/xlsx twins)."""
    matrices = []
    companions = {}
    for proj in sorted(os.listdir(DATA)):
        pdir = project_data_dir(proj)
        if not os.path.isdir(pdir) or not proj.startswith("RPMM"):
            continue
        cands = []
        for root, _dirs, files in os.walk(pdir):
            for fn in sorted(files):
                p = os.path.join(root, fn)
                rel = logical_path(p)
                if SKIP_RE.search(fn):
                    continue
                if COMPANION_RE.search(fn) or re.search(r"alkane", fn, re.I):
                    # alkane series = RI calibration standard runs, not
                    # per-metabolite quantification -> provenance companion
                    companions.setdefault(proj, []).append(rel)
                    continue
                if re.search(r"pheno", fn, re.I):
                    companions.setdefault(proj, []).append(rel)
                    continue
                if not re.search(r"\.(txt|csv|tsv|xlsx?|xls)$", fn, re.I):
                    continue
                head = ""
                if not re.search(r"\.xlsx?$|\.xls$", fn, re.I):
                    try:
                        with open(p, encoding="utf-8", errors="replace") as f:
                            head = f.readline().strip()
                    except OSError:
                        continue
                    if head.startswith("Name:"):          # MSP or alkane series
                        companions.setdefault(proj, []).append(rel)
                        continue
                cands.append(p)
        # de-duplicate: when foo.txt and foo.xlsx coexist, keep the text one
        by_stem = {}
        for p in cands:
            stem = os.path.splitext(os.path.basename(p))[0].lower()
            by_stem.setdefault(stem, []).append(p)
        for stem, group in by_stem.items():
            text = [p for p in group if not re.search(r"\.xlsx?$|\.xls$", p, re.I)]
            chosen = group[0] if not text else min(text, key=len)
            matrices.append(chosen)
            for p in group:
                if p != chosen:
                    companions.setdefault(proj, []).append(
                        logical_path(p))
    return matrices, companions


def project_titles():
    import json
    p = PROJECTS_JSON
    if not os.path.exists(p):
        return {}
    with open(p) as f:
        items = json.load(f)
    return {i["p"].rsplit("/", 1)[-1]: i.get("label", "") for i in items}


def cat_columns(pcols):
    """Phenodata columns usable as study variables (categorical, compact)."""
    out = []
    for name, vals in pcols.items():
        v = [x for x in vals if x and x.lower() not in ("null", "na")]
        u = sorted(set(v))
        if 2 <= len(u) <= 8 and all(len(x) <= 40 for x in u) and \
                not re.search(r"species|binomial|weight|order", name, re.I):
            out.append(name)
    return out


def pick_companion(rels, rel_src, predicate):
    """Same-directory companion wins; else the project-unique match."""
    same = [r for r in rels
            if predicate(r) and os.path.dirname(r) == os.path.dirname(rel_src)]
    if same:
        return same[0]
    anyp = [r for r in rels if predicate(r)]
    return anyp[0] if anyp else ""


def project_pheno_values(proj):
    """All values appearing in the project's phenodata files (sample-ID pool)."""
    vals = set()
    for root, _dirs, files in os.walk(project_data_dir(proj)):
        for fn in sorted(files):
            if re.search(r"pheno", fn, re.I) and not re.search(r"\.xlsx?$|\.xls$", fn, re.I):
                try:
                    with open(os.path.join(root, fn), encoding="utf-8",
                              errors="replace") as f:
                        for ln in f.read().splitlines():
                            for cell in ln.replace(",", "\t").split("\t"):
                                cell = cell.strip().strip('"')
                                if cell:
                                    vals.add(cell)
                except OSError:
                    pass
    return vals


SML_DESTS = {"chemical_name", "feature_id", "unique_peak_id", "window_token",
             "retention_index", "mass_to_charge"}


def looks_like_sample_metadata(mx, pheno_vals):
    """A table whose rows describe samples (e.g. a run-order workbook), not
    per-metabolite quantification: no SML annotation columns at all, and its
    meta name/sample column values match the project's sample IDs."""
    produced = set()
    for key in mx["meta"]:
        for src, dst in SML_META_MAP:
            if key == norm_name(src):
                produced.add(dst)
                break
    if produced & SML_DESTS:
        return False
    for key in ("sample_name", "name", "sample", "sample_id", "biological_name",
                "gcid"):
        vals = mx["meta"].get(key)
        if vals and pheno_vals and \
                sum(1 for v in vals if v in pheno_vals) >= max(3, 0.5 * len(vals)):
            return True
    # preparation-tracking sheets (run order, tube weights, solvent volumes):
    # sample-name column, zero metabolite annotation, and every "sample"
    # column is a weight/volume/order descriptor rather than a sample ID
    name_col = next((mx["meta"].get(k) for k in
                     ("sample_name", "name", "sample", "sample_id",
                      "biological_name", "gcid") if mx["meta"].get(k)), None)
    if name_col and mx["samples"] and all(
            re.search(r"weight|volume|order|combined|抽出|分析|加量", s, re.I)
            for s in mx["samples"]):
        return True
    return False


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

# matrix annotation columns -> SML column names (first match wins)
SML_META_MAP = [
    ("numbering_for_matrix", "feature_id"),
    ("numbering_for_unique_peaks", "unique_peak_id"),
    ("metabolites", "chemical_name"),
    ("final_annotation", "chemical_name"),
    ("synonym", "chemical_name"),
    ("annotation", "chemical_name"),
    ("compound", "chemical_name"),
    ("name", "chemical_name"),
    ("iis_i_a_or_u", "annotation_grade"),
    ("i_a_m_or_p", "annotation_type"),
    ("ri", "retention_index"),
    ("win", "window_token"),
    ("ret_umin", "retention_time"),
    ("rt", "retention_time"),
    ("rt_min", "retention_time"),
    ("retention_time", "retention_time"),
    ("mass", "mass_to_charge"),
    ("mz", "mass_to_charge"),
    ("m_z", "mass_to_charge"),
]


def build_sml_table(mx):
    """Annotation columns in SML naming + abundance_<sample> columns."""
    n = mx["n_rows"]
    cols = {"sml_index": list(range(1, n + 1))}
    used = {}
    for src, dst in SML_META_MAP:
        key = norm_name(src)
        if key in mx["meta"] and dst not in used:
            vals = mx["meta"][key]
            if dst in ("retention_index", "retention_time", "mass_to_charge"):
                cols[dst] = [to_num(v) for v in vals]
            elif dst == "window_token":
                raw = vals
                toks = [re.sub(r"^.*?(Win[0-9]+_\w+).*$", r"\1", v)
                        if re.search(r"Win[0-9]+_", v) else (v or None)
                        for v in raw]
                cols["window_token"] = toks
                # RT in seconds lives as the leading number of the Win cell
                rt = []
                for v in raw:
                    m = re.match(r"\s*([0-9.]+)", v or "")
                    rt.append(float(m.group(1)) if m else math.nan)
                cols.setdefault("retention_time_seconds", rt)
            else:
                cols[dst] = [v or None for v in vals]
            used[key] = dst
    for key, vals in mx["meta"].items():
        if key in used:
            continue
        dst = "meta_" + key
        numeric = sum(1 for v in vals[:60] if not math.isnan(to_num(v)))
        cols[dst] = ([to_num(v) for v in vals]
                     if numeric >= 0.5 * min(60, len(vals)) else
                     [v or None for v in vals])
    for s in mx["samples"]:
        cols["abundance_" + s] = mx["abund"][s]
    return cols


def write_parquet(path, cols):
    arrays, names = [], []
    for name, vals in cols.items():
        nonnull = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
        if not nonnull or all(isinstance(v, int) for v in nonnull):
            t, data = pa.int64(), [None if v is None or (isinstance(v, float) and math.isnan(v)) else int(v) for v in vals]
        elif all(isinstance(v, (int, float)) for v in nonnull):
            t, data = pa.float64(), [None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v) for v in vals]
        else:
            t, data = pa.string(), [None if v is None else str(v) for v in vals]
        arrays.append(pa.array(data, type=t))
        names.append(name)
    pq.write_table(pa.Table.from_arrays(arrays, names=names), path)


def export_matrix(mx, title, species_raw, species_cv, species_src,
                  pheno_rel, msp_rel, out_dir, curate):
    stem = os.path.splitext(os.path.basename(mx["path"]))[0]
    proj = logical_path(mx["path"]).split(os.sep)[0]
    base = f"{proj}_{stem}"
    rel_src = logical_path(mx["path"])
    # provenance rows: every pair carries its source files into the log
    for kind, rel in (("source_matrix", rel_src),
                      ("source_phenodata", pheno_rel),
                      ("source_msp", msp_rel)):
        if rel:
            curate.append((proj, base, kind, rel))
    pheno_path = physical_path(pheno_rel) if pheno_rel else None
    pcols = read_phenodata(pheno_path) if pheno_path else None
    id_col, aligned = align_phenodata(pcols, mx["samples"]) if pcols else (None, {})
    if pcols and len(aligned) < 0.5 * len(mx["samples"]):
        curate.append((proj, base, "phenodata_partial_match",
                       f"{len(aligned)}/{len(mx['samples'])} samples matched "
                       f"via '{id_col}'"))
    unmatched = [s for s in mx["samples"] if s not in aligned]
    if unmatched and len(unmatched) <= 12:
        curate.append((proj, base, "samples_without_phenodata",
                       ", ".join(unmatched)))
    elif unmatched:
        curate.append((proj, base, "samples_without_phenodata",
                       f"{len(unmatched)} samples: " +
                       ", ".join(unmatched[:8]) + ", ..."))
    for name in mx["unknown_columns"]:
        curate.append((proj, base, "unclassified_column_kept_as_meta",
                       f"'{name}' is text and unrecognized; exported as "
                       f"meta_{norm_name(name)} (not a sample column)"))

    groups = cat_columns(pcols) if pcols else []
    polarity = ("negative" if re.search(r"neg|negative", stem, re.I)
                else "positive")

    # ---- TOML: sample metadata for the mzTab-M MTD table ----
    loc_base = ("http://metabobank.riken.jp/data/" +
                os.path.dirname(rel_src).replace(os.sep, "/"))
    samples_tbl = []
    for i, s in enumerate(mx["samples"], start=1):
        e = {
            "index": i,
            "name": s,
            "assay": s,
            "ms_run_ref": f"ms_run[{i}]",
            "location": f"{loc_base}/{os.path.basename(rel_src)}#{s}",
            "description": f"{proj} sample {s}",
        }
        if aligned.get(s):
            for k, v in sorted(aligned[s].items()):
                if v:
                    e[norm_name(k) or f"col{k}"] = v
        samples_tbl.append(e)

    toml = {
        "schema_version": 1,
        "generated_by": "investigation/export_tables.py",
        "generated_on": GENERATED_ON,
        "mzd_table": "MTD",
        "project": {
            "id": proj,
            "title": title or proj,
            "external_study": f"http://metabobank.riken.jp/pmm/{proj}",
        },
        "matrix": {
            "name": stem,
            "source_file": rel_src,
            "format": mx["format"],
            "n_features": mx["n_rows"],
            "n_samples": len(mx["samples"]),
            "metadata_columns": mx["meta_columns"],
            "unclassified_columns": mx["unknown_columns"],
        },
        "companion_files": {
            "phenodata": logical_path(pheno_path) if pheno_path else "",
            "msp": msp_rel or "",
        },
        "organism": {
            "raw": species_raw or "",
            "cv": ("[NCBITaxon, %s, %s, ]" % species_cv if species_cv else ""),
            "source": species_src or "",
        },
        "instrument": {
            "technique": "GC-MS" if "win" in mx["meta"] else "LC-MS",
            "scan_polarity": polarity,
            "scan_polarity_source": "derived from filename; confirm before MTD use",
        },
        "sample": samples_tbl,
    }
    # unique study-variable values observed
    uvals = {}
    for g in groups:
        uvals[g] = sorted({e.get(norm_name(g)) for e in samples_tbl
                           if e.get(norm_name(g))})
    toml["study_variables"] = {"candidate_columns": groups,
                               "levels": uvals}

    toml_path = os.path.join(out_dir, base + ".toml")
    with open(toml_path, "wb") as f:
        f.write(tomli_w.dumps(toml).encode())

    # ---- Parquet: per-metabolite quantification (SML) ----
    parq_path = os.path.join(out_dir, base + ".parquet")
    write_parquet(parq_path, build_sml_table(mx))

    return {
        "project": proj, "matrix": stem, "source": rel_src,
        "toml": output_reference(toml_path), "parquet": output_reference(parq_path),
        "n_features": mx["n_rows"], "n_samples": len(mx["samples"]),
        "phenodata": logical_path(pheno_path) if pheno_path else "",
        "id_col": id_col or "", "samples_matched": len(aligned),
        "species_cv": species_cv[1] if species_cv else "",
        "study_variables": ", ".join(groups),
    }


def main():
    global DATA, OUT, REPORTS, GENERATED_ON
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("projects", nargs="*", help="optional RPMM project IDs")
    parser.add_argument("--projects-root", type=Path, default=REPO / "projects",
                        help="root containing RPMM####/data input directories")
    parser.add_argument("--output-root", type=Path,
                        help="artifact root containing RPMM####/tables (default: projects root)")
    parser.add_argument("--reports-root", type=Path,
                        help="global report directory (default: reports/tables; explicit for subsets)")
    parser.add_argument("--generated-on", help="ISO date for reproducible TOML metadata")
    args = parser.parse_args()
    DATA = args.projects_root.resolve()
    OUT = (args.output_root or DATA).resolve()
    REPORTS = (args.reports_root or REPO / "reports" / "tables").resolve()
    if args.projects and args.reports_root is None:
        parser.error("Project subsets require --reports-root so they cannot replace full reports")
    date = args.generated_on
    if date is None and os.environ.get("SOURCE_DATE_EPOCH"):
        date = datetime.datetime.fromtimestamp(int(os.environ["SOURCE_DATE_EPOCH"]),
                                               datetime.timezone.utc).date().isoformat()
    if date is not None:
        GENERATED_ON = datetime.date.fromisoformat(date).isoformat()
    only = args.projects
    os.makedirs(REPORTS, exist_ok=True)
    titles = project_titles()
    matrices, companions = inventory()
    if only:
        matrices = [m for m in matrices
                    if logical_path(m).split(os.sep)[0] in only]
    print(f"inventory: {len(matrices)} matrices over "
          f"{len({logical_path(m).split(os.sep)[0] for m in matrices})} projects")

    curate, report = [], []
    errors = []
    for m in sorted(matrices):
        proj = logical_path(m).split(os.sep)[0]
        try:
            mx = read_matrix(m)
            rel_src = logical_path(m)
            rels = companions.get(proj, [])
            pheno_rel = pick_companion(
                rels, rel_src,
                lambda r: bool(re.search(r"pheno", os.path.basename(r), re.I)))
            msp_rel = pick_companion(
                rels, rel_src,
                lambda r: bool(re.search(r"nist", os.path.basename(r), re.I))
                and r.lower().endswith(".txt"))
            try:
                pcols_probe = (read_phenodata(physical_path(pheno_rel))
                               if pheno_rel else None)
            except Exception:
                pcols_probe = None
            # guard: sample-metadata tables (run-order workbooks etc.) are not
            # per-metabolite quantification -> companion, not an export
            if looks_like_sample_metadata(mx, project_pheno_values(proj)):
                companions.setdefault(proj, []).append(logical_path(m))
                curate.append((proj, os.path.basename(m),
                               "sample_metadata_table_skipped",
                               "rows describe samples (match project sample IDs), "
                               "no metabolite/feature annotation; exported nowhere"))
                print(f"  skip {logical_path(m)}  (sample metadata table)")
                continue
            raw, src = find_species_raw(os.path.dirname(m), proj, pcols_probe)
            cv = resolve_species(raw) if raw else None
            if raw and not cv:
                curate.append((proj, os.path.basename(m), "species_unresolved_taxon",
                               f"raw '{raw}' has no curated NCBITaxon mapping; "
                               "organism.cv left empty"))
            if not raw:
                curate.append((proj, os.path.basename(m), "species_missing",
                               "no species string in MSP or phenodata"))
            out_dir = os.path.join(OUT, proj, "tables")
            os.makedirs(out_dir, exist_ok=True)
            rep = export_matrix(mx, titles.get(proj, ""), raw,
                                cv, src, pheno_rel, msp_rel, out_dir, curate)
            report.append(rep)
            print(f"  ok  {rep['source']}  ({rep['n_features']} features x "
                  f"{rep['n_samples']} samples)")
        except Exception as e:
            errors.append((proj, logical_path(m), str(e)))
            curate.append((proj, os.path.basename(m), "conversion_failed", str(e)))
            print(f"  ERR {logical_path(m)}: {e}")

    rep_path = os.path.join(REPORTS, "_report.tsv")
    fields = ["project", "matrix", "source", "toml", "parquet", "n_features",
              "n_samples", "phenodata", "id_col", "samples_matched",
              "species_cv", "study_variables"]
    with open(rep_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(report)
    cur_path = os.path.join(REPORTS, "_curation.tsv")
    with open(cur_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["project", "subject", "issue", "detail"])
        w.writerows(curate)
    print(f"\nreport:   {rep_path}  ({len(report)} matrices)")
    print(f"curation: {cur_path}  ({len(curate)} issues)")
    if errors:
        print(f"FAILED:   {len(errors)} matrices")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
