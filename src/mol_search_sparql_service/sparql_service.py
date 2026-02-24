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
def similarity_search(
    smiles: str,
    limit: int = 10,
    db_names: str | None = None,
    fp_type: str = "morgan_ecfp",
    use_chirality: bool = False,
    min_score: float = 0.0,
) -> list[SearchResult]:
    """Perform similarity search using precomputed fingerprints.

    Args:
        smiles: Query SMILES string.
        limit: Maximum number of results to return.
        db_names: Optional database name to filter results.
        fp_type: Fingerprint type key to use.
        use_chirality: Whether to respect chirality when computing fingerprints.
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
            use_chirality=use_chirality,
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
    smart: str,
    limit: int = 100,
    db_names: str | None = None,
    use_chirality: bool = False,
    min_match_count: int = 1,
) -> list[SubstructureSearchResult]:
    """Perform substructure search using a SMARTS/SMILES pattern.

    Args:
        smart: Query SMARTS or SMILES pattern to match.
        limit: Maximum number of results to return (default: 100).
        db_names: Optional database name to limit the search.
        use_chirality: Whether to respect chirality for matching.
        min_match_count: Minimum number of substructure matches required.

    Example:
        ```sparql
        PREFIX func: <urn:sparql-function:>
        SELECT ?result ?matchCount WHERE {
            [] a func:SubstructureSearch ;
                func:smart "c1ccccc1" ;
                func:result ?result ;
                func:matchCount ?matchCount .
        }
        ```
    """
    try:
        db_list = [db_names] if db_names else None
        results = engine.search_substructure(
            smart,
            limit=limit,
            db_names=db_list,
            use_chirality=use_chirality,
            min_match_count=min_match_count,
        )
        return [
            SubstructureSearchResult(result=URIRef(r.id), matchCount=int(r.match_count))
            for r in results
        ]
    except Exception as e:
        print(f"Error in substructure_search: {e}")
        return []


def generate_docs() -> str:
    """Return markdown documentation for all SPARQL functions and fingerprint options."""
    # g.generate_docs() returns markdown for every registered @g.type_function
    docs = g.generate_docs()
    # Append fingerprint options table
    fp_lines = ["## Available Fingerprint Types\n"]
    fp_lines.append("| Key | Short Name | Description |")
    fp_lines.append("|-----|------------|-------------|")
    for key, fp in FINGERPRINTS.items():
        desc = " ".join(fp.description.split())
        fp_lines.append(f"| `{key}` | {fp.short_name} | {desc} |")
    return "\n## Functions\n\n" + docs.rstrip() + "\n\n" + "\n".join(fp_lines) + "\n"


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
    return f"""You are an expert in writing SPARQL queries for the Chemistry Search Service.
Use the following documentation to write correct queries that strictly match the defined `func:` namespace and classes:

{generate_docs()}

The SPARQL endpoint accepts `HTTP GET` and `HTTP POST`, typically at `/sparql` on the hosting server.
"""


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup: load the compounds file into the engine if needed (multi-worker case)
    compounds_file = os.environ.get("COMPOUNDS_FILE")
    if compounds_file and not engine.datasets:
        print(f"Initializing engine from: {compounds_file}")
        try:
            engine.load_file(compounds_file)
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
