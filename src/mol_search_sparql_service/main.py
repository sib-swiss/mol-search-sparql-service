import argparse
import uvicorn
import os
import sys
import socket

from .rdkit_fingerprints import engine, FINGERPRINTS
from .sparql_service import app


def _validate_fingerprint_types(fp_types: list[str] | None) -> None:
    """Validate that all specified fingerprint types exist. Dies with error message if not."""
    if not fp_types:
        return
    invalid_types = [fp for fp in fp_types if fp not in FINGERPRINTS]
    if invalid_types:
        available = ", ".join(sorted(FINGERPRINTS.keys()))
        parser.error(
            f"Unknown fingerprint type(s): {', '.join(invalid_types)}\n"
            f"Available types: {available}"
        )


def _validate_port(port: int) -> None:
    """Validate port is in valid range and available. Dies with error message if not."""
    if not (1 <= port <= 65535):
        parser.error(f"Port must be between 1 and 65535, got {port}")

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(("localhost", port)) == 0:
            parser.error(
                f"Port {port} is already in use. Use -p/--port to specify a different port."
            )


def main() -> None:
    global parser
    parser = argparse.ArgumentParser(description="Start the Chemistry SPARQL Service")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--file", help="Path to compounds.tsv")
    group.add_argument("-s", "--sparql", help="Path to SPARQL query file")

    parser.add_argument(
        "-e", "--endpoint", help="SPARQL endpoint URL (required if using -s)"
    )
    parser.add_argument(
        "-p",
        "--port",
        type=int,
        default=8010,
        help="Port to run the server on (default: 8010)",
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help="Number of Uvicorn workers (default: 1)",
    )
    parser.add_argument(
        "-d",
        "--daemon",
        action="store_true",
        help="Daemonize the server process (after fingerprints have been computed and all validation checks passed)",
    )
    parser.add_argument(
        "-t",
        "--fingerprints",
        type=str,
        help="Comma-separated list of fingerprint types to compute (e.g. morgan_ecfp,pattern). If omitted, all types are computed.",
    )
    args = parser.parse_args()

    # === VALIDATION PHASE (before any irreversible operations) ===
    if args.sparql and not args.endpoint:
        parser.error("-s/--sparql requires -e/--endpoint")

    if args.file and not os.path.exists(args.file):
        parser.error(f"Compounds file not found: {args.file}")

    if args.sparql and not os.path.exists(args.sparql):
        parser.error(f"SPARQL query file not found: {args.sparql}")

    fp_types = args.fingerprints.split(",") if args.fingerprints else None
    _validate_fingerprint_types(fp_types)

    if args.fingerprints:
        os.environ["FINGERPRINTS_LIST"] = args.fingerprints

    _validate_port(args.port)

    # === DATA LOADING PHASE (errors still visible to user) ===
    # 1. Load Data
    if args.file:
        engine.load_file(args.file, fp_types=fp_types)
    else:
        with open(args.sparql, "r") as f:
            query = f.read()
        engine.load_from_sparql(args.endpoint, query, fp_types=fp_types)

    # === DAEMONIZATION PHASE (only after all checks and preparations complete) ===
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
        with open(os.devnull, "r") as f:
            os.dup2(f.fileno(), sys.stdin.fileno())
        with open("server.log", "a") as f:
            os.dup2(f.fileno(), sys.stdout.fileno())
            os.dup2(f.fileno(), sys.stderr.fileno())

    # 2. Start Server
    print(
        f"Starting SPARQL endpoint on port {args.port} with {args.workers} worker(s)..."
    )
    import rdflib.plugins.sparql

    print(
        f"DEBUG: CUSTOM_EVALS keys: {list(rdflib.plugins.sparql.CUSTOM_EVALS.keys())}"
    )

    # Uvicorn requires an import string when using multiple workers
    if args.workers > 1:
        uvicorn.run(
            "mol_search_sparql_service.sparql_service:app",
            host="0.0.0.0",
            port=args.port,
            workers=args.workers,
        )
    else:
        uvicorn.run(app, host="0.0.0.0", port=args.port)


if __name__ == "__main__":
    main()
