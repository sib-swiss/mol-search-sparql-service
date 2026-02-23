from dataclasses import dataclass
from typing import List, Any
from rdflib import Literal, URIRef, Namespace
from rdflib_endpoint import SparqlEndpoint, DatasetExt
from mcp.server.fastmcp import FastMCP

from .rdkit_fingerprints import engine, FINGERPRINTS

# Define Namespace
FUNC = Namespace("urn:sparql-function:")

# Define Result Dataclass
@dataclass
class SubstructureSearchResult:
    result: URIRef
    matchCount: int

@dataclass
class SearchResult:
    result: URIRef
    score: float

@dataclass
class FingerprintInfo:
    fpType: str
    description: str
    mechanism: str
    shortName: str

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

@g.type_function(FUNC)
def list_fingerprints() -> List[FingerprintInfo]:
    """
### func:ListFingerprints
Lists available fingerprint types.
- `func:fpType` (string, output): The fingerprint type identifier.
- `func:description` (string, output): Description of the fingerprint.
- `func:shortName` (string, output): Short name (e.g., ECFP).
- `func:mechanism` (string, output): Explanation of how it works.
    """
    return [
        FingerprintInfo(
            fpType=key,
            description=val.get('description', ''),
            mechanism=val.get('explainability', {}).get('mechanism', ''),
            shortName=val.get('short_name', '')
        )
        for key, val in FINGERPRINTS.items()
    ]

@g.type_function(FUNC)
def similarity_search(smiles: str, limit: int = 10, db_names: str = None, fp_type: str = 'morgan_ecfp', use_chirality: bool = False, min_score: float = 0.0) -> List[SearchResult]:
    """
### func:SimilaritySearch
Performs similarity search based on fingerprints.
- `func:smiles` (string, required): Query SMILES string.
- `func:limit` (integer, optional): Maximum results (default 10).
- `func:dbNames` (string, optional): Filter by database source.
- `func:fpType` (string, optional): Fingerprint type (default 'morgan_ecfp').
- `func:useChirality` (boolean, optional): Whether to respect chirality (default false).
- `func:minScore` (float, optional): Minimum similarity score (default 0.0).
- `func:result` (URI, output): The matching compound URI.
- `func:score` (float, output): Tanimoto similarity score (0-1).
    """
    try:
        if fp_type not in FINGERPRINTS:
            print(f"Error: Invalid fingerprint type '{fp_type}'")
            return []
            
        db_list = [db_names] if db_names else None
        results = engine.search_similarity(smiles, limit=limit, db_names=db_list, fp_type=fp_type, use_chirality=use_chirality, min_score=min_score)
        return [SearchResult(result=URIRef(r['compound']['id']), score=float(r['similarity'])) for r in results]
    except Exception as e:
        print(f"Error in similarity_search: {e}")
        return []

@g.type_function(FUNC)
def substructure_search(smart: str, limit: int = 100, db_names: str = None, use_chirality: bool = False, min_match_count: int = 1) -> List[SubstructureSearchResult]:
    """
### func:SubstructureSearch
Performs substructure search.
- `func:smart` (string, required): Query SMARTS or SMILES pattern.
- `func:limit` (integer, optional): Maximum results (default 100).
- `func:dbNames` (string, optional): Filter by database source.
- `func:useChirality` (boolean, optional): Whether to respect chirality (default false).
- `func:minMatchCount` (integer, optional): Minimum matches required (default 1).
- `func:result` (URI, output): The matching compound URI.
- `func:matchCount` (integer, output): Number of matches found (1 if boolean match).
    """
    try:
        db_list = [db_names] if db_names else None
        results = engine.search_substructure(smart, limit=limit, db_names=db_list, use_chirality=use_chirality, min_match_count=min_match_count)
        return [SubstructureSearchResult(result=URIRef(r['id']), matchCount=int(r.get('match_count', 1))) for r in results]
    except Exception as e:
        print(f"Error in substructure_search: {e}")
        return []


# =========================================================================
# MCP Documenting & Helper Endpoints
# =========================================================================
# The MCP server does not provide functional search tools bypassing SPARQL.
# Instead, it provides documentation resources and prompts so that LLMs
# can easily understand how to construct valid SPARQL against the service.
# =========================================================================

mcp = FastMCP("Chemistry Search")

# Dynamic Docstring Injection
def _update_docstrings():
    fp_descriptions = []
    for key, val in FINGERPRINTS.items():
        desc = f"- {key}: {val['description']}"
        fp_descriptions.append(desc)
    
    fp_doc = "\n    ".join(fp_descriptions)
    
    if similarity_search.__doc__:
        similarity_search.__doc__ += f"\n    Available Fingerprint Types:\n    {fp_doc}\n    "

_update_docstrings()

def _assemble_schema() -> str:
    header = """
# SPARQL Schema for Chemistry Search Service

## Namespace
Prefix: `func:` <urn:sparql-function:>

## Classes
"""
    classes = [
        similarity_search.__doc__,
        substructure_search.__doc__,
        list_fingerprints.__doc__
    ]
    footer = """
## Examples

### Identify Similar Molecules
```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result ?score WHERE {
    [] a func:SimilaritySearch ;
       func:smiles "c1ccccc1" ;
       func:fpType "morgan_ecfp" ;
       func:result ?result ;
       func:score ?score .
}
```

### Identify Substructures
```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result WHERE {
    [] a func:SubstructureSearch ;
       func:smart "C(=O)O" ;
       func:result ?result .
}
```
"""
    return header + "\n".join(d.strip() for d in classes if d) + footer

@mcp.resource("sparql://schema")
def sparql_schema() -> str:
    """Returns the generated SPARQL schema for the chemistry search service."""
    return _assemble_schema()

@mcp.prompt()
def sparql_assistant() -> str:
    """
    Interactive prompt providing the SPARQL schema and writing guidelines for the Chemistry Search Service.
    """
    return f"""You are an expert in writing SPARQL queries for the specific Chemistry Search Service.
Use the following documentation to assist in writing correct queries strictly matching the defined `func:` namespaces and classes:

{_assemble_schema()}

The SPARQL endpoint listens to `HTTP GET` and `HTTP POST` typically located at `/sparql` on the hosting server.
"""


# Create the SparqlEndpoint using the DatasetExt 'g'
app = SparqlEndpoint(
    graph=g,
    path="/sparql",
    cors_enabled=True,
    # Functions are already registered via decorators on 'g'
)

# Mount MCP Server (FastMCP's internal app)
app.mount("/mcp", mcp.sse_app())
