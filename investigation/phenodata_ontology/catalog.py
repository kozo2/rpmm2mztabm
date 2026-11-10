"""Verified, offline ontology vocabulary used by the Phenodata annotator.

Normal annotation runs only read the checked-in JSON snapshot.  ``--refresh``
explicitly fetches official sources, resolves GitHub files to immutable commits,
and records source hashes.  OLS responses themselves are retained as evidence
because the live API is not an immutable release endpoint.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import hashlib
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request


DATA = Path(__file__).resolve().parent / "data"
OBO_REPOS = {
    "PO": ("Planteome/plant-ontology", "po.obo"),
    "PATO": ("pato-ontology/pato", "pato.obo"),
    "UO": ("bio-ontology-research-group/unit-ontology", "uo.obo"),
    "PECO": ("Planteome/plant-experimental-conditions-ontology", "peco.obo"),
}
API_TERMS = {
    "obi": ["OBI:0100026"],
    "cob": ["COB:0000022"],
    "efo": ["EFO:0000635", "EFO:0000513", "EFO:0000727", "EFO:0000399",
            "EFO:0000724", "EFO:0000688", "EFO:0000324", "EFO:0000428"],
    "chebi": ["CHEBI:15366", "CHEBI:15377", "CHEBI:26710", "CHEBI:46793",
              "CHEBI:41879", "CHEBI:749470"],
    "ncbitaxon": ["NCBITaxon:3702", "NCBITaxon:4081", "NCBITaxon:4530",
                  "NCBITaxon:4236", "NCBITaxon:3983", "NCBITaxon:3827",
                  "NCBITaxon:3847", "NCBITaxon:3708", "NCBITaxon:3711",
                  "NCBITaxon:4678", "NCBITaxon:35875", "NCBITaxon:130426",
                  "NCBITaxon:4683", "NCBITaxon:47664", "NCBITaxon:1402",
                  "NCBITaxon:1390", "NCBITaxon:210367"],
}


def canonical_id(identifier: str) -> str:
    """Convert a standard OBO IRI/short form to a CURIE without guessing IDs."""
    identifier = identifier.strip()
    if "/obo/" in identifier:
        identifier = identifier.rsplit("/", 1)[-1]
    elif "/efo/" in identifier and identifier.rsplit("/", 1)[-1].startswith("EFO_"):
        identifier = identifier.rsplit("/", 1)[-1]
    return re.sub(r"^([A-Za-z]+)_(\d+)$", r"\1:\2", identifier)


class OntologyCatalog:
    """Read-only, versioned vocabulary. Unknown IDs fail closed."""

    def __init__(self, data_dir: Path | str | None = None):
        self.data_dir = Path(data_dir) if data_dir else DATA
        term_bytes = (self.data_dir / "ontology_terms.json").read_bytes()
        term_data = json.loads(term_bytes)
        source_data = json.loads((self.data_dir / "ontology_sources.json").read_text())
        if _sha(term_bytes) != source_data["terms_sha256"]:
            raise ValueError("Ontology snapshot checksum does not match ontology_sources.json")
        for source in source_data["sources"].values():
            if source.get("response_utf8") is not None:
                if _sha(source["response_utf8"].encode("utf-8")) != source["sha256"]:
                    raise ValueError(f"OLS response checksum does not match: {source['url']}")
        self.terms = term_data["terms"]
        self.sources = source_data["sources"]
        self.snapshot_created_at = source_data["created_at"]

    def get(self, identifier: str) -> dict | None:
        return self.terms.get(canonical_id(identifier))

    def is_current(self, identifier: str) -> bool:
        term = self.get(identifier)
        return bool(term and term.get("verified") and not term["obsolete"]
                    and term.get("label") and term.get("source_id") in self.sources)

    def validate(self, identifier: str) -> dict:
        """Return an explicit validation result; never replace deprecated IDs."""
        term = self.get(identifier)
        return {
            "id": canonical_id(identifier),
            "valid": self.is_current(identifier),
            "reason": ("unknown_id" if term is None else
                       "obsolete_id" if term["obsolete"] else "verified"),
            "replaced_by": term.get("replaced_by", []) if term else [],
        }

    def find_exact(self, text: str, ontology: str | None = None) -> list[dict]:
        """Return current label/EXACT-synonym candidates, not final assignments."""
        text = " ".join(text.strip().lower().split())
        matches = []
        for term in self.terms.values():
            if ontology and term["ontology"].lower() != ontology.lower():
                continue
            if not self.is_current(term["id"]):
                continue
            labels = [term["label"]] + [s["label"] for s in term["synonyms"]
                                        if s["scope"] == "EXACT"]
            if text in {" ".join(label.lower().split()) for label in labels}:
                matches.append(term)
        return matches


_DEFAULT_CATALOG: OntologyCatalog | None = None


def _default() -> OntologyCatalog:
    global _DEFAULT_CATALOG
    if _DEFAULT_CATALOG is None:
        _DEFAULT_CATALOG = OntologyCatalog()
    return _DEFAULT_CATALOG


def get(identifier: str) -> dict | None:
    """Module-level convenience wrapper for the checked-in snapshot."""
    return _default().get(identifier)


def is_current(identifier: str) -> bool:
    return _default().is_current(identifier)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _download(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "rpmm-explorer-phenodata/1.0"})
    with urllib.request.urlopen(req, timeout=45) as response:
        return response.read()


def _quoted(line: str) -> str:
    match = re.search(r'"((?:[^"\\]|\\.)*)"', line)
    return re.sub(r"\\(.)", r"\1", match[1]) if match else ""


def _synonym_scope(scope: str) -> str:
    scope = scope.rsplit("#", 1)[-1].upper()
    for value in ("EXACT", "BROAD", "NARROW", "RELATED"):
        if scope in (value, "HAS" + value + "SYNONYM"):
            return value
    return "RELATED"


def _parse_obo(data: bytes, prefix: str, source_id: str, source: dict) -> dict:
    text = data.decode("utf-8")
    terms = {}
    for stanza in re.split(r"^\[Term\]\s*$", text, flags=re.MULTILINE)[1:]:
        stanza = re.split(r"^\[", stanza, maxsplit=1, flags=re.MULTILINE)[0]
        fields = {}
        for line in stanza.splitlines():
            if ": " in line:
                key, value = line.split(": ", 1)
                fields.setdefault(key, []).append(value)
        term_id = fields.get("id", [""])[0]
        if not term_id.startswith(prefix + ":") or not fields.get("name"):
            continue
        synonyms = []
        for value in fields.get("synonym", []):
            match = re.search(r'"\s+(EXACT|BROAD|NARROW|RELATED)\b', value)
            synonyms.append({"label": _quoted(value), "scope": match[1] if match else "RELATED"})
        terms[term_id] = {
            "id": term_id, "label": fields["name"][0],
            "iri": "http://purl.obolibrary.org/obo/" + term_id.replace(":", "_"),
            "ontology": prefix,
            "definition": _quoted(fields.get("def", [""])[0]),
            "obsolete": fields.get("is_obsolete") == ["true"],
            "replaced_by": [canonical_id(v.split()[0]) for v in fields.get("replaced_by", [])],
            "consider": [canonical_id(v.split()[0]) for v in fields.get("consider", [])],
            "alt_ids": fields.get("alt_id", []), "synonyms": synonyms,
            "parents": [v.split()[0] for v in fields.get("is_a", [])],
            "relationships": fields.get("relationship", []),
            "source_id": source_id, "source_url": source["resolved_url"],
            "ontology_version": source["version"], "verified": True,
        }
    return terms


def _fetch_obo(item: tuple[str, tuple[str, str]], timestamp: str) -> tuple[dict, dict]:
    prefix, (repo, filename) = item
    commit_url = (f"https://api.github.com/repos/{repo}/commits?"
                  + urllib.parse.urlencode({"path": filename, "per_page": 1}))
    commit_response = _download(commit_url)
    commit = json.loads(commit_response)[0]["sha"]
    url = f"https://raw.githubusercontent.com/{repo}/{commit}/{filename}"
    data = _download(url)
    version_match = re.search(rb"^data-version: (.+)$", data, re.MULTILINE)
    source_id = f"{prefix}@{commit}"
    source = {
        "ontology": prefix, "kind": "official_obo", "resolved_url": url,
        "url": f"https://raw.githubusercontent.com/{repo}/master/{filename}",
        "version": version_match[1].decode() if version_match else commit,
        "commit": commit, "retrieved_at": timestamp, "sha256": _sha(data),
        "bytes": len(data), "commit_metadata_url": commit_url,
        "commit_metadata_sha256": _sha(commit_response),
    }
    return _parse_obo(data, prefix, source_id, source), {source_id: source}


def _fetch_api_term(ontology: str, identifier: str, info: dict, timestamp: str) -> tuple[dict, dict]:
    url = (f"https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}/terms?"
           + urllib.parse.urlencode({"obo_id": identifier}))
    data = _download(url)
    records = json.loads(data).get("_embedded", {}).get("terms", [])
    records = [r for r in records if r.get("obo_id") == identifier]
    if len(records) != 1:
        raise ValueError(f"Expected exactly one OLS record for {identifier}: {len(records)}")
    record = records[0]
    source_id = "OLS:" + identifier
    prefix = identifier.split(":")[0]
    config = info.get("config", {})
    source = {
        "ontology": prefix, "kind": "official_ols_api_snapshot", "url": url,
        "resolved_url": url, "version": config.get("version") or config.get("versionIri", ""),
        "version_iri": config.get("versionIri", ""), "retrieved_at": timestamp,
        "sha256": _sha(data), "bytes": len(data),
        "response_utf8": data.decode("utf-8"),
        "ontology_metadata_url": f"https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}",
    }
    synonyms = record.get("obo_synonym") or []
    # OLS's flattened synonyms omit synonym scope; retain those as RELATED.
    parsed_synonyms = [{"label": s.get("name", ""), "scope": _synonym_scope(str(s.get("scope", "RELATED")))}
                       for s in synonyms if isinstance(s, dict)]
    included = {s["label"] for s in parsed_synonyms}
    parsed_synonyms += [{"label": s, "scope": "RELATED"}
                        for s in record.get("synonyms", []) if s not in included]
    replacement = record.get("term_replaced_by")
    if not replacement:
        replacement = record.get("annotation", {}).get("term replaced by", [])
    if isinstance(replacement, str):
        replacement = [replacement]
    term = {
        "id": identifier, "label": record["label"], "iri": record["iri"],
        "ontology": prefix, "definition": "\n".join(record.get("description") or []),
        "obsolete": bool(record.get("is_obsolete")),
        "replaced_by": [canonical_id(x) for x in replacement or []],
        "consider": [], "alt_ids": [], "synonyms": parsed_synonyms,
        "parents": [], "relationships": [], "source_id": source_id,
        "source_url": url, "ontology_version": source["version"], "verified": True,
    }
    return {identifier: term}, {source_id: source}


def refresh(additional_ids: list[str] | None = None) -> dict:
    """Fetch a new explicit snapshot. Fail on unavailable sources or invalid IDs."""
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
    terms, sources = {}, {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        for new_terms, new_sources in pool.map(lambda i: _fetch_obo(i, timestamp), OBO_REPOS.items()):
            terms.update(new_terms)
            sources.update(new_sources)
    requests = {k: list(v) for k, v in API_TERMS.items()}
    for identifier in additional_ids or []:
        prefix = identifier.split(":")[0].lower()
        if prefix in {p.lower() for p in OBO_REPOS}:
            if identifier not in terms:
                raise ValueError(f"ID absent from official OBO: {identifier}")
            continue
        requests.setdefault(prefix, [])
        if identifier not in requests[prefix]:
            requests[prefix].append(identifier)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        infos = dict(zip(requests, pool.map(
            lambda ontology: json.loads(_download(f"https://www.ebi.ac.uk/ols4/api/ontologies/{ontology}")),
            requests)))
        jobs = [(ontology, identifier) for ontology, identifiers in requests.items() for identifier in identifiers]
        for new_terms, new_sources in pool.map(
                lambda job: _fetch_api_term(job[0], job[1], infos[job[0]], timestamp), jobs):
            terms.update(new_terms)
            sources.update(new_sources)
    DATA.mkdir(parents=True, exist_ok=True)
    term_document = {"schema_version": 1, "terms": dict(sorted(terms.items()))}
    encoded = (json.dumps(term_document, ensure_ascii=False, indent=2) + "\n").encode()
    (DATA / "ontology_terms.json").write_bytes(encoded)
    source_document = {"schema_version": 1, "created_at": timestamp,
                       "terms_sha256": _sha(encoded), "sources": dict(sorted(sources.items()))}
    (DATA / "ontology_sources.json").write_text(
        json.dumps(source_document, ensure_ascii=False, indent=2) + "\n")
    return {"terms": len(terms), "sources": len(sources), "terms_sha256": _sha(encoded)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="Explicitly retrieve a new official vocabulary snapshot")
    parser.add_argument("--add-term", action="append", default=[])
    args = parser.parse_args()
    if args.refresh:
        print(json.dumps(refresh(args.add_term), indent=2))
    else:
        catalog = OntologyCatalog()
        print(json.dumps({"terms": len(catalog.terms), "sources": len(catalog.sources),
                          "created_at": catalog.snapshot_created_at}, indent=2))
