import os
from dataclasses import dataclass
from contextlib import asynccontextmanager
from rdflib import URIRef, Namespace
from rdflib_endpoint import SparqlEndpoint, DatasetExt
from fastapi import FastAPI
from mcp.server.fastmcp import FastMCP

from .rdkit_fingerprints import engine, FINGERPRINTS

# Define Namespace
FUNC = Namespace("urn:sparql-function:")


# Define Result Dataclass
@dataclass
class SubstructureSearchResult:
    result: URIRef
    """The URI of the matching compound."""
    matchCount: int
    """Number of matches found (1 if boolean match)."""
    matchedSmiles: str
    """SMILES of the matched fragment, rendered from the target (stereo preserved)."""
    matchedSmarts: str
    """SMARTS of the matched fragment, rendered from the target (stereo preserved)."""


@dataclass
class SearchResult:
    result: URIRef
    """The URI of the matching compound."""
    score: float
    """Tanimoto similarity score (0-1)."""


@dataclass
class FingerprintInfo:
    fpType: str
    """Identifier key for the fingerprint type (e.g. `morgan_ecfp`)."""
    description: str
    """Human readable description of the fingerprint."""
    mechanism: str
    """Explainability mechanism / how bits map to substructures."""
    shortName: str
    """Short display name for the fingerprint (e.g. ECFP, MACCS)."""


@dataclass
class DatabaseInfo:
    dbName: str
    """The name or URI of the database stored in the service."""


# Initialize DatasetExt
g = DatasetExt()
g.bind("func", FUNC)


# =========================================================================
# SPARQL ENDPOINTS (rdflib_endpoint format)
# =========================================================================
# The following bindings exist specifically to satisfy the `SparqlEndpoint`
# library interfaces, expecting strict typings and returning specific class
# dataclasses that are serialized over RDF logic.
# =========================================================================


@g.type_function()
def list_fingerprints() -> list[FingerprintInfo]:
    """List available fingerprint types.

    Example:
        ```sparql
        PREFIX func: <urn:sparql-function:>
        SELECT ?fpType ?description ?shortName WHERE {
            [] a func:ListFingerprints ;
                func:fpType ?fpType ;
                func:description ?description ;
                func:shortName ?shortName .
        }
        ```
    """
    return [
        FingerprintInfo(
            fpType=key,
            description=val.description,
            mechanism=val.explainability.mechanism,
            shortName=val.short_name,
        )
        for key, val in FINGERPRINTS.items()
    ]


@g.type_function()
def list_databases() -> list[DatabaseInfo]:
    """List available database names loaded in the service.

    Example:
        ```sparql
        PREFIX func: <urn:sparql-function:>
        SELECT ?dbName WHERE {
            [] a func:ListDatabases ;
                func:dbName ?dbName .
        }
        ```
    """
    return [DatabaseInfo(dbName=db) for db in engine.get_databases()]


@g.type_function()
def similarity_search(
    smiles: str,
    limit: int = 10,
    db_names: str | None = None,
    fp_type: str = "morgan_ecfp",
    min_score: float = 0.0,
) -> list[SearchResult]:
    """Perform similarity search using precomputed fingerprints.

    Args:
        smiles: Query SMILES string.
        limit: Maximum number of results to return.
        db_names: Optional database name to filter results.
        fp_type: Fingerprint type key to use.
        min_score: Minimum similarity score threshold (0.0 - 1.0).

    Example:
        ```sparql
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?score WHERE {
            [] a func:SimilaritySearch ;
                func:smiles "[NH3+][C@@H](Cc1ccccc1)C(=O)[O-]" ;
                func:limit 3 ;
                func:result ?result ;
                func:score ?score .
        }
        ```
    """
    try:
        if fp_type not in FINGERPRINTS:
            print(f"Error: Invalid fingerprint type '{fp_type}'")
            return []

        db_list = [db_names] if db_names else None
        results = engine.search_similarity(
            smiles,
            limit=limit,
            db_names=db_list,
            fp_type=fp_type,
            min_score=min_score,
        )
        return [
            SearchResult(result=URIRef(r.compound.id), score=float(r.similarity))
            for r in results
        ]
    except Exception as e:
        print(f"Error in similarity_search: {e}")
        return []


@g.type_function()
def substructure_search(
    smart: str | None = None,
    smiles: str | None = None,
    limit: int = 100,
    db_names: str | None = None,
    min_match_count: int = 1,
    use_chirality: bool = False,
) -> list[SubstructureSearchResult]:
    """Perform substructure search using a SMARTS or SMILES query pattern.

    Provide exactly one of ``func:smart`` or ``func:smiles``; the two are parsed
    differently by RDKit:

    - ``func:smart`` is parsed with ``MolFromSmarts`` as a query graph. It keeps
      query features (wildcards, any-bond ``~``, recursive SMARTS, degree/charge
      constraints) and matches literally — an aliphatic ``C`` will NOT match an
      aromatic carbon, and aromaticity is not re-perceived.
    - ``func:smiles`` is parsed with ``MolFromSmiles``: sanitized, aromatized and
      with stereochemistry perceived, then used as a substructure query.

    Args:
        smart: Query SMARTS pattern to match (mutually exclusive with smiles).
        smiles: Query SMILES pattern to match (mutually exclusive with smart).
        limit: Maximum number of results to return (default: 100).
        db_names: Optional database name to limit the search.
        min_match_count: Minimum number of substructure matches required.
        use_chirality: If true, both tetrahedral (R/S) and double-bond (E/Z) stereochemistry are enforced during matching. Defaults to false. Most meaningful for SMILES queries; SMARTS encodes its own stereo in the pattern.

    Example (SMARTS):
        ```sparql
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?matchCount ?matchedSmiles ?matchedSmarts WHERE {
            [] a func:SubstructureSearch ;
                func:smart "[#6]~[#7]" ;
                func:result ?result ;
                func:matchCount ?matchCount ;
                func:matchedSmiles ?matchedSmiles ;
                func:matchedSmarts ?matchedSmarts .
        }
        ```

    Example (SMILES):
        ```sparql
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?matchCount ?matchedSmiles ?matchedSmarts WHERE {
            [] a func:SubstructureSearch ;
                func:smiles "c1ccccc1" ;
                func:result ?result ;
                func:matchCount ?matchCount ;
                func:matchedSmiles ?matchedSmiles ;
                func:matchedSmarts ?matchedSmarts .
        }
        ```
    """
    try:
        if smart and smiles:
            print("Error: provide either func:smart or func:smiles, not both")
            return []
        if smart:
            query, query_type = smart, "smarts"
        elif smiles:
            query, query_type = smiles, "smiles"
        else:
            print("Error: substructure_search requires func:smart or func:smiles")
            return []

        db_list = [db_names] if db_names else None
        results = engine.search_substructure(
            query,
            limit=limit,
            db_names=db_list,
            min_match_count=min_match_count,
            use_chirality=use_chirality,
            query_type=query_type,
        )
        # Emit one row per distinct matched fragment so each match's SMILES and
        # SMARTS are individually bindable. Compounds with no renderable fragment
        # still surface once with empty fragment strings.
        rows = []
        for r in results:
            fragments = list(zip(r.matched_smiles, r.matched_smarts)) or [("", "")]
            for smi, sma in fragments:
                rows.append(
                    SubstructureSearchResult(
                        result=URIRef(r.id),
                        matchCount=int(r.match_count),
                        matchedSmiles=smi,
                        matchedSmarts=sma,
                    )
                )
        return rows
    except Exception as e:
        print(f"Error in substructure_search: {e}")
        return []


def generate_docs() -> str:
    """Return markdown documentation for all SPARQL functions and fingerprint options."""
    # g.generate_docs() returns markdown for every registered @g.type_function
    docs = g.generate_docs()
    # Append fingerprint options table
    fp_lines = ["## 🫆 Available Fingerprint Types\n"]
    fp_lines.append("| Key | Short Name | Description |")
    fp_lines.append("|-----|------------|-------------|")
    for key, fp in FINGERPRINTS.items():
        desc = " ".join(fp.description.split())
        fp_lines.append(f"| `{key}` | {fp.short_name} | {desc} |")
    return "\n## 📖 Functions\n\n" + docs.rstrip() + "\n\n" + "\n".join(fp_lines) + "\n"


# =========================================================================
# MCP Documenting & Helper Endpoints
# =========================================================================
# The MCP server does not provide functional search tools bypassing SPARQL.
# Instead, it provides documentation resources and prompts so that LLMs
# can easily understand how to construct valid SPARQL against the service.
# =========================================================================

mcp = FastMCP("Chemistry Search SPARQL service")


@mcp.resource("sparql://schema")
def sparql_schema() -> str:
    """Returns the generated SPARQL schema for the chemistry search service."""
    return generate_docs()


@mcp.prompt()
def sparql_assistant() -> str:
    """Interactive prompt providing schema and guidelines for writing SPARQL against this service."""
    sparql_url = os.environ.get("SPARQL_PUBLIC_URL", "/sparql")
    return f"""You are an expert in writing SPARQL queries for the Chemistry Search Service.
Use the following documentation to write correct queries that strictly match the defined `func:` namespace and classes:

{generate_docs()}

The SPARQL endpoint accepts `HTTP GET` and `HTTP POST` at: {sparql_url}
"""


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup: load the compounds file into the engine if needed (multi-worker case)
    compounds_file = os.environ.get("COMPOUNDS_FILE")
    if compounds_file and not engine.datasets:
        print(f"Initializing engine from: {compounds_file}")

        fp_list = os.environ.get("FINGERPRINTS_LIST")
        fp_types = fp_list.split(",") if fp_list else None

        # Empty/unset CACHE_DIR means caching is disabled
        cache_dir = os.environ.get("CACHE_DIR") or None

        try:
            engine.load_file(compounds_file, fp_types=fp_types, cache_dir=cache_dir)
        except Exception as e:
            print(f"Error loading compounds file on startup: {e}")
    yield


app = SparqlEndpoint(
    graph=g,
    path="/sparql",
    cors_enabled=True,
    lifespan=lifespan,
)


# Mount MCP Server (FastMCP's internal app)
app.mount("/mcp", mcp.sse_app())
