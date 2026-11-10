"""Input manifest provenance and optional, cached PMM RDF discovery.

Discovery records real file links only. It never manufactures sample URIs from
CSV IDs or propagates a project-level species to all measurements.
"""

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

DATA = Path(__file__).parent / "data"
ENDPOINT = "https://metabobank.riken.jp/pmm/sparql"
GRAPH = "http://metadb.riken.jp/db/plantMetabolomics"


def manifest_records(investigation):
    path = investigation / "manifest.csv"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    return {row["local_path"].removeprefix("data/") if hasattr(str, "removeprefix") else
            row["local_path"][5:]: row for row in records}


def discovery_query(file_iris):
    iris = " ".join("<" + iri + ">" for iri in sorted(set(file_iris)))
    return ("SELECT ?file ?parent ?predicate WHERE { GRAPH <" + GRAPH +
            "> { VALUES ?file { " + iris + " } ?parent ?predicate ?file . } }")


def refresh_discovery(file_iris, timeout=20):
    """Persist successful response OR a truthful retrieval-failure record."""
    query = discovery_query(file_iris)
    query_sha = hashlib.sha256(query.encode()).hexdigest()
    request_url = ENDPOINT + "?" + urllib.parse.urlencode({
        "query": query, "format": "application/sparql-results+json"})
    record = {"endpoint": ENDPOINT, "graph": GRAPH, "query": query,
              "query_sha256": query_sha, "retrieved_at_utc": datetime.now(timezone.utc).isoformat()}
    try:
        request = urllib.request.Request(request_url, headers={
            "Accept": "application/sparql-results+json",
            "User-Agent": "rpmm-explorer-phenodata-annotation/1.0"})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            raw = response.read()
        result = json.loads(raw)
        if "results" not in result:
            raise ValueError("Endpoint did not return SPARQL SELECT results")
        record.update(status="retrieved", response_sha256=hashlib.sha256(raw).hexdigest(),
                      response=result,
                      limitation="Parent-link discovery only; sample-level joins must be verified before annotation reuse")
    except Exception as error:
        record.update(status="unavailable", error=f"{type(error).__name__}: {error}",
                      limitation="No sample-level RDF annotations were reused; source CSV evidence is used instead")
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "rdf_discovery.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    return record


def load_discovery():
    path = DATA / "rdf_discovery.json"
    if path.exists():
        return json.loads(path.read_text())
    return {"status": "not_requested", "endpoint": ENDPOINT,
            "limitation": "No verified sample-level RDF context cache is available"}
