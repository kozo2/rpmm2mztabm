"""Download every File downloadURL in downloadurls_final4.tsv.

- follows redirects
- verifies the MD5 recorded in the endpoint against the bytes received
- stores each remote project under projects/RPMM####/data/
- preserves data/RPMM####/... as the logical manifest path
- writes manifest.csv summarising every file
"""

import argparse
import concurrent.futures as cf
import csv
import hashlib
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = (HERE.parent / "projects").resolve()

# The endpoint records "AraMetLeaves.Rdata" but the server only serves ".RData".
URL_FIXES = {
    "http://metabobank.riken.jp/data/RPMM0001/AnalysedDataset/AraMetLeaves.Rdata":
        "http://metabobank.riken.jp/data/RPMM0001/AnalysedDataset/AraMetLeaves.RData",
}


def logical_path(url):
    """Keep the original remote data/RPMM####/... identity in the manifest."""
    path = urllib.parse.unquote(urllib.parse.urlparse(url).path)
    path = path.rstrip()                      # drops the stray trailing space
    return path.lstrip("/")


def local_path(url):
    """Map data/RPMM####/... to projects/RPMM####/data/... safely."""
    parts = Path(logical_path(url)).parts
    if (len(parts) < 3 or parts[0] != "data" or not parts[1].startswith("RPMM")
            or not parts[1][4:].isdigit()):
        raise ValueError(f"unsupported RPMM download path: {url}")
    target = (ROOT / parts[1] / "data" / Path(*parts[2:])).resolve()
    if not str(target).startswith(str(ROOT) + os.sep):
        raise ValueError(f"refusing path outside download root: {target}")
    target.relative_to((ROOT / parts[1] / "data").resolve())
    return target


def fetch(row):
    url = row["downloadURL"]
    real = URL_FIXES.get(url, url)
    dest = local_path(url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out = {
        **row, "requested": real, "local_path": logical_path(url),
        "status": "", "bytes": 0, "md5_expected": row.get("md5", ""),
        "md5_actual": "", "md5_match": "", "error": "",
    }
    if dest.exists() and dest.stat().st_size > 0:
        pass  # re-verify below
    try:
        req = urllib.request.Request(real, headers={"User-Agent": "curl/8.7"})
        with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as fh:
            h = hashlib.md5()
            n = 0
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
                h.update(chunk)
                n += len(chunk)
        out["status"] = "ok"
        out["bytes"] = n
        out["md5_actual"] = h.hexdigest()
        exp = row.get("md5", "")
        out["md5_match"] = ("yes" if exp == out["md5_actual"] else
                            "no" if exp else "no-md5-in-endpoint")
    except urllib.error.HTTPError as e:
        out["status"] = f"HTTP {e.code}"
        out["error"] = str(e)
    except Exception as e:  # noqa: BLE001
        out["status"] = "error"
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def main():
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projects-root", type=Path, default=HERE.parent / "projects")
    parser.add_argument("--urls", type=Path, default=HERE / "downloadurls_final4.tsv")
    parser.add_argument("--manifest", type=Path, default=HERE / "manifest.csv")
    args = parser.parse_args()
    ROOT = args.projects_root.resolve()
    with args.urls.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        parser.error("download URL input is empty")
    print(f"rows: {len(rows)}  ->  {ROOT}", flush=True)
    results = []
    with cf.ThreadPoolExecutor(max_workers=6) as ex:
        for i, res in enumerate(ex.map(fetch, rows), 1):
            results.append(res)
            if i % 25 == 0 or i == len(rows):
                ok = sum(1 for r in results if r["status"] == "ok")
                print(f"  [{i}/{len(rows)}] ok={ok}", flush=True)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
        w.writeheader()
        w.writerows(results)

    ok = [r for r in results if r["status"] == "ok"]
    bad = [r for r in results if r["status"] != "ok"]
    match = [r for r in ok if r["md5_match"] == "yes"]
    mismatch = [r for r in ok if r["md5_match"] == "no"]
    print(f"\nDONE  ok={len(ok)}  failed={len(bad)}")
    print(f"bytes: {sum(r['bytes'] for r in ok):,}")
    print(f"md5 matched: {len(match)}  mismatched: {len(mismatch)}")
    for r in bad:
        print(f"  FAIL {r['status']} {r['downloadURL']} {r['error'][:80]}")
    for r in mismatch:
        print(f"  MISMATCH {r['downloadURL']} expected={r['md5_expected']} got={r['md5_actual']}")


if __name__ == "__main__":
    main()
