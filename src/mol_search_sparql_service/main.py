import argparse
import uvicorn
import requests
import tempfile
import os
import sys

from .rdkit_fingerprints import engine
from .sparql_service import app

def main():
    parser = argparse.ArgumentParser(description="Start the Chemistry SPARQL Service")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--file', help='Path to compounds.tsv')
    group.add_argument('-s', '--sparql', help='Path to SPARQL query file')

    parser.add_argument('-e', '--endpoint', help='SPARQL endpoint URL (required if using -s)')
    parser.add_argument('-p', '--port', type=int, default=8010, help='Port to run the server on (default: 8010)')
    parser.add_argument('-w', '--workers', type=int, default=1, help='Number of Uvicorn workers (default: 1)')
    parser.add_argument('-d', '--daemon', action='store_true', help='Daemonize the server process')
    args = parser.parse_args()

    if args.sparql and not args.endpoint:
        parser.error("-s/--sparql requires -e/--endpoint")

    # 1. Compile Data
    if args.file:
        engine.load_file(args.file)
    else:
        with open(args.sparql, 'r') as f:
            query = f.read()

        print(f"Fetching data from {args.endpoint}...")
        try:
            # We request TSV directly from the SPARQL endpoint
            headers = {'Accept': 'text/tab-separated-values'}
            resp = requests.post(args.endpoint, data={'query': query}, headers=headers)
            resp.raise_for_status()

            fd, temp_path = tempfile.mkstemp(suffix='.tsv')
            # Important: Keep the temp file around for workers by not deleting it immediately
            with os.fdopen(fd, 'wb') as f:
                f.write(resp.content)

            # Pass to environment so workers can use it
            os.environ['COMPOUNDS_FILE'] = temp_path
            os.environ['DELETE_COMPOUNDS_FILE'] = '1'

            engine.load_file(temp_path)
            # We will rely on sparql_service or OS to clean this up, or clean it up after uvicorn exits.
        except Exception as e:
            print(f"Failed to fetch data from endpoint: {e}")
            sys.exit(1)

    # Store explicit file path in environment for workers
    if args.file:
        os.environ['COMPOUNDS_FILE'] = args.file

    # Daemonize if requested
    if args.daemon:
        # Double fork to detach from terminal completely
        if os.fork() > 0:
            sys.exit(0)
        os.setsid()
        if os.fork() > 0:
            sys.exit(0)

        # Redirect standard file descriptors
        sys.stdout.flush()
        sys.stderr.flush()
        with open(os.devnull, 'r') as f:
            os.dup2(f.fileno(), sys.stdin.fileno())
        with open('server.log', 'a') as f:
            os.dup2(f.fileno(), sys.stdout.fileno())
            os.dup2(f.fileno(), sys.stderr.fileno())

    # 2. Start Server
    print(f"Starting SPARQL endpoint on port {args.port} with {args.workers} worker(s)...")
    import rdflib.plugins.sparql
    print(f"DEBUG: CUSTOM_EVALS keys: {list(rdflib.plugins.sparql.CUSTOM_EVALS.keys())}")

    # Uvicorn requires an import string when using multiple workers
    if args.workers > 1:
        uvicorn.run("mol_search_sparql_service.sparql_service:app", host="0.0.0.0", port=args.port, workers=args.workers)
    else:
        uvicorn.run(app, host="0.0.0.0", port=args.port)

if __name__ == "__main__":
    main()
