from typing import Any

import pytest
import requests
import time
import socket
from multiprocessing import Process, set_start_method
import uvicorn
from mol_search_sparql_service.rdkit_fingerprints import engine
from mol_search_sparql_service.sparql_service import app

PORT = 8000
URL = f"http://localhost:{PORT}/sparql"

def is_port_in_use(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def _run_server():
    print("Loading data for test SPARQL server in child process...")
    engine.load_file("compounds.tsv")
    uvicorn.run(app, host="0.0.0.0", port=PORT, log_level="info")

@pytest.fixture(scope="module", autouse=True)
def sparql_server():
    server_process = None
    if not is_port_in_use(PORT):
        print(f"\nStarting test SPARQL server on port {PORT}...")
        server_process = Process(
            target=_run_server,
            daemon=True,
        )
        server_process.start()
        # Wait for server to start
        start_time = time.time()
        timeout = 20

        while time.time() - start_time < timeout:
            if not server_process.is_alive():
                raise RuntimeError("Server process crashed during startup.")

            try:
                # Basic check to see if we get a response
                requests.get(f"http://localhost:{PORT}/sparql?query=SELECT * WHERE {{}} LIMIT 1", timeout=1)
                break
            except requests.exceptions.RequestException:
                time.sleep(1)
        else:
            server_process.kill()
            server_process.join()
            raise RuntimeError("Server failed to start in time.")
    yield
    if server_process:
        print("\nStopping test SPARQL server...")
        server_process.kill()
        server_process.join()


def sparql_query(query: str) -> Any:
    response = requests.get(URL, params={'query': query}, headers={'Accept': 'application/sparql-results+json'})
    response.raise_for_status()
    return response.json()["results"]["bindings"]



# --- TESTS ---

def test_similarity_search_default():
    test_mol = "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]"
    bindings = sparql_query(f"""
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?score WHERE {{
            [] a func:SimilaritySearch ;
                func:smiles "{test_mol}" ;
                func:result ?result ;
                func:score ?score .
        }}
    """)
    assert len(bindings) > 0
    assert "score" in bindings[0]


def test_similarity_search_limit():
    test_mol = "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]"
    bindings = sparql_query(f"""
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?score WHERE {{
            [] a func:SimilaritySearch ;
                func:smiles "{test_mol}" ;
                func:limit 3 ;
                func:result ?result ;
                func:score ?score .
        }}
    """)
    assert len(bindings) > 0
    assert len(bindings) <= 3


def test_substructure_search_default():
    # Benzene ring
    bindings = sparql_query(f"""
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?matchCount WHERE {{
            [] a func:SubstructureSearch ;
                func:smart "c1ccccc1" ;
                func:result ?result ;
                func:matchCount ?matchCount .
        }}
    """)
    assert len(bindings) > 0
    assert "matchCount" in bindings[0]

def test_list_fingerprints():
    bindings = sparql_query("""
        PREFIX func: <urn:sparql-function:>
        SELECT ?fpType ?description ?shortName WHERE {
            [] a func:ListFingerprints ;
                func:fpType ?fpType ;
                func:description ?description ;
                func:shortName ?shortName .
        }
        """)
    assert len(bindings) > 0

    fp_types = [b["fpType"]["value"] for b in bindings]
    assert "morgan_ecfp" in fp_types
    assert "maccs" in fp_types
