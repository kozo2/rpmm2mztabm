# PMM File Download Report

**Source list:** [`downloadurls_final4.txt`](downloadurls_final4.txt) — 299 URLs (after excluding `RawDataset`, `DropMet`, `/ABF/`, `/MTBLS45/`, `/StandardMixtures/`, `EachBatch_raw.zip`)
**Current destination:** [`../projects/`](../projects/) — `RPMM####/data/` 内に元の内部ディレクトリ構造を保持。取得時の論理パス `data/RPMM####/...` はmanifest内で維持しています。
**Manifest:** [`manifest.csv`](manifest.csv)
**Date:** 2026-10-15

## Result

| Metric | Value |
|---|---:|
| URLs requested | 299 |
| Files downloaded | **299 / 299** |
| Total size | **120,772,741 bytes (115.2 MiB)** |
| MD5 verified against endpoint | **288** |
| MD5 mismatch (stale endpoint metadata) | **11** |
| Failed downloads | **0** |

All 10 `.xlsx` archives pass `zipfile.testzip()`; all 6 `.xls` files carry the OLE2 magic
`d0cf11e0a1b11ae1`; the single `.rdata` is correctly gzip-compressed (`1f8b08`).

## Two malformed URLs in the endpoint (worked around)

1. `…/RPMM0024/AnalysedDataset/MetMat_abs.txt%20` → **404**.
   The trailing `%20` is spurious; the real file is `MetMat_abs.txt` (**200**, 333,759 bytes).
   Once corrected, its MD5 matches the endpoint exactly (`eabc56fca76d0d3f10322f581b7226ef`),
   confirming only the *stored URL* is wrong, not the data.

2. `…/RPMM0001/AnalysedDataset/AraMetLeaves.Rdata` → **404** (case-sensitive filesystem).
   `AraMetLeaves.RData` serves fine (9,130 bytes), but its bytes do **not** match the recorded
   MD5 — see below.

## 11 MD5 mismatches

`curl` and `urllib` fetch byte-identical content for every one of these, so this is **not** a
download error: the endpoint's recorded MD5s are stale or cross-wired.

| URL | Note |
|---|---|
| `RPMM0018/Year2/gc012/AnalysedDataset/pheno.csv` | **MD5s swapped with the next row** |
| `RPMM0018/Year2/gc2/AnalysedDataset/pheno.csv` | **MD5s swapped with the previous row** |
| `RPMM0001/AnalysedDataset/AraMetLeaves.Rdata` | only served as `.RData` |
| `RPMM0001/AnalysedDataset/12918_2007_53_MOESM6_ESM.xls` | content valid OLE2 |
| `RPMM0014/AnalysedDataset/GCMS/UV-integrated.1day.xls` | content valid OLE2 |
| `RPMM0014/AnalysedDataset/GCMS/UV-integrated.4day.xls` | content valid OLE2 |
| `RPMM0020/AnalysedDataset/CRMN_normalized.txt` | content valid |
| `RPMM0023/AnalysedDataset/080522arabi_seeds_nv_AnnotCombMat.txt` | content valid |
| `RPMM0037/AnalysedDataset/phenodata.csv` | content valid |
| `RPMM0041/AnalysedDataset/080814Miyaosan_sample02_nv_MetMat_marked.xls` | content valid OLE2 |
| `RPMM0068/AnalysedDataset/LC/Phenodata.csv` | content valid |

The RPMM0018 pair is a genuine **metadata cross-wiring**: the MD5 the endpoint records for
`gc012/pheno.csv` (`e190fc…`) is the MD5 of the `gc2` file actually served, and vice versa.
The two files differ (14,503 vs 14,399 bytes), so both downloads are real distinct files —
only the endpoint's hash assignments are transposed.

## Reproducing

```bash
python3 investigation/download_all.py     # download + verify (6 threads)
```

`download_all.py` folds in the `.Rdata` → `.RData` fix via its `URL_FIXES` map and validates
the destination path stays inside the selected `projects/RPMM####/data/` folder before writing.

上記の件数と検証結果は取得当時の記録です。再配置後の実ファイル対応とハッシュは
[`../reports/relocation/manifest.json`](../reports/relocation/manifest.json) に保存しています。
