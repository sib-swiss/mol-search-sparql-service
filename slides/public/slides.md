## Extending SPARQL with Python

**Make Python functions queryable in a SPARQL endpoint**

Marco Pagni, Vincent Emonet · RDF Focus Group

---

## SPARQL is great for graphs

SPARQL excels at traversal, pattern matching, and federated queries across RDF knowledge graphs.

But knowledge graphs are not just data, they are increasingly paired with computation:

- Cheminformatics: molecular similarity, fingerprinting, substructure search
- NLP: named entity recognition, embedding similarity
- ML inference: scoring, classification, ranking

Python is one of the most popular solution for in research and data libraries


> **SPARQL has no native way to run arbitrary computation.**

[github.com/w3c/sparql-dev/issues/131](https://github.com/w3c/sparql-dev/issues/131)

---

## What are our options today?

**1. SPARQL extension functions**: most triplestores support `BIND(ext:myFunc(...) AS ?x)`, but registration is store-specific, requires deep internals knowledge, and admin control over the endpoint.

**2. Post-processing**: fetch all results, compute in Python, lose the ability to filter or join server-side. Pulls too much data.

**3. Federated `SERVICE` endpoint**: seems the right architecture, but building a SPARQL-compliant HTTP server from scratch is a lot of plumbing.

---

## Inspirations

- Goes back 10 years ago with Jerven's identifiers.org custom java service
- FastAPI for the API design and dev experience
- RDFLib for RDF in python

Extending an existing known library: `rdflib-endpoint` (by Vincent)

RDFLib is a RDF python library with everything needed to do most RDF task (e.g. SPARQL query engine)

RDFLib-endpoint enables to easily deploy graph built with RDFLib

Already used by projects like [bioregistry.io/sparql](https://bioregistry.io/sparql)

---

## The goal

Write a Python function. Get a working SPARQL endpoint. No triplestore plugins, no RDF parsing, no boilerplate.

```python
from rdflib_endpoint import DatasetExt

ds = DatasetExt()

@dataclass
class SearchResult:
    result: URIRef
    score: float

@ds.type_function()
def similarity_search(smiles: str, limit: int = 10) -> list[SearchResult]:
    results = engine.search(smiles, limit=limit)
    return [SearchResult(result=URIRef(r.id), score=r.score) for r in results]
```

Queryable immediately:

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?compound ?score WHERE {
    [] a func:SimilaritySearch ;
       func:smiles "c1ccccc1" ;
       func:result ?compound ;
       func:score ?score .
}
```

---

## `DatasetExt`: 4 decorator patterns

`rdflib-endpoint` extends RDFLib's `Dataset` with decorator helpers, covering most custom evaluation use-cases observed in the wild:

| Decorator             | Triggered by SPARQL pattern          | Best for                        |
| --------------------- | ------------------------------------ | ------------------------------- |
| `@type_function`      | `[] a func:FuncName ; func:arg ?v`   | Multi-input / multi-output      |
| `@predicate_function` | `?s func:predName ?o`                | Object computed from subject    |
| `@extension_function` | `BIND(func:funcName(...) AS ?var)`   | Scalar / multi-binding results  |
| `@graph_function`     | `BIND(func:funcName(...) AS ?g)`     | Return a temporary named graph  |

Python snake_case names are automatically mapped to SPARQL camelCase/PascalCase IRIs under the configured namespace (default to `urn:sparql-function:`).

---

## How rdflib-endpoint works

`rdflib-endpoint` builds on RDFLib's custom evaluation hook to intercept query patterns and dispatch them to Python functions.

**The flow:**

1. Create a RDFLib `ds = DatasetExt()`
2. Decorate a Python function with `@ds.type_function()` (or other patterns)
3. rdflib-endpoint registers it in RDFLib's evaluation engine
4. Serve the dataset as a SPARQL 1.1 endpoint (using FastAPI)
5. Incoming queries that match the pattern trigger your function
6. Return values become SPARQL result bindings automatically

Idiomatic python: type annotations and dataclasses handle all the IRI and binding mapping.

---

## `@type_function`: rich triple pattern

The most expressive pattern. Inputs and outputs are predicates in the triple pattern.

```python
@dataclass
class SearchResult:
    result: URIRef
    """The URI of the matching compound."""
    score: float
    """Tanimoto similarity score (0-1)."""

@ds.type_function()
def similarity_search(smiles: str, limit: int = 10) -> list[SearchResult]:
    """Similarity search over precomputed fingerprints.
    
    Args:
        smiles: Query SMILES string.
        limit: Maximum number of results to return.

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
        \```
    """
    return [SearchResult(result=URIRef(r.id), score=r.score)
            for r in engine.search(smiles, limit=limit)]
```

Each field in the dataclass becomes a SPARQL predicate automatically.

Proper docstring help generating better docs and populate Yasgui tabs automatically.

----



```python
@dataclass
class SearchResult:
    result: URIRef
    """The URI of the matching compound."""
    score: float
    """Tanimoto similarity score (0-1)."""

@ds.type_function()
def similarity_search(smiles: str, limit: int = 10) -> list[SearchResult]:
    """Similarity search over precomputed fingerprints.
    
    Args:
        smiles: Query SMILES string.
        limit: Maximum number of results to return.

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
    return [SearchResult(result=URIRef(r.id), score=r.score)
            for r in engine.search(smiles, limit=limit)]
```

---

## `@predicate_function`: looks like a triple

The function fires when its IRI appears as a predicate. The subject is the input, the return value is the object.

```python
import bioregistry
from rdflib import OWL, URIRef

conv = bioregistry.get_converter()

@ds.predicate_function(namespace=OWL._NS)
def same_as(input_iri: URIRef) -> list[URIRef]:
    """Return all equivalent IRIs via Bioregistry."""
    prefix, identifier = conv.compress(input_iri).split(":", 1)
    return [URIRef(iri) for iri in bioregistry.get_providers(prefix, identifier).values()]
```

```sparql
SELECT ?sameAs WHERE {
    <https://identifiers.org/CHEBI/1> owl:sameAs ?sameAs .
}
```

Looks like ordinary triple resolution -- no `BIND`, no `SERVICE` wrapper needed.

---

## `@extension_function`: classic SPARQL extension

The familiar `BIND(func:name(...) AS ?var)` pattern. Return a list to emit multiple rows, use a dataclass to populate multiple variables.

```python
@dataclass
class SplitResult:
    value: str
    index: int

@ds.extension_function()
def split_index(input_str: str, separator: str = ",") -> list[SplitResult]:
    return [SplitResult(value=p, index=i) for i, p in enumerate(input_str.split(separator))]
```

```sparql
SELECT ?input ?part ?partIndex WHERE {
    VALUES ?input { "hello world" "cheese is good" }
    BIND(func:splitIndex(?input, " ") AS ?part)
}
```

`?part` gets `SplitResult.value`, `?partIndex` gets `SplitResult.index` -- variable names derived from dataclass field names.

> Additional variables work like labels on Wikidata label service: concatenate to the bound `?variable`

---

## `@graph_function`: populate a temporary graph

New special concept!

---

## Serving the endpoint

One line to deploy as a standalone app:

```python
from rdflib_endpoint import SparqlEndpoint

app = SparqlEndpoint(graph=ds, title="My SPARQL Service")
```

Run:

```sh
uvicorn main:app --reload
```

Or mount as a router inside an existing FastAPI/Flask/Django app:

```python
app.include_router(SparqlRouter(graph=ds, path="/sparql"))
```

SPARQL query examples embedded in docstrings are surfaced automatically as YASGUI tabs.

---

## Implementation: chemistry search

[`mol-search-sparql-service`](https://github.com/sib-swiss/mol-search-sparql-service): RDKit fingerprint search exposed entirely via `@type_function`

| SPARQL type               | What it does                               |
| ------------------------- | ------------------------------------------ |
| `func:ListFingerprints`   | Enumerate available fingerprint algorithms |
| `func:SimilaritySearch`   | Tanimoto similarity over preloaded dataset |
| `func:SubstructureSearch` | SMARTS/SMILES substructure match           |

---

## Federation: call it from anywhere

Because the endpoint speaks standard SPARQL 1.1, any other endpoint can delegate to it with `SERVICE`.

```sparql
#+ endpoint: https://sparql.rhea-db.org/sparql
PREFIX func: <urn:sparql-function:>
SELECT ?compound ?score WHERE {
    ?compound a <http://rdf.rhea-db.org/Compound> .
    SERVICE <https://biosoda.unil.ch/mol-search/sparql> {
        [] a func:SimilaritySearch ;
           func:smiles "c1ccccc1" ;
           func:result ?compound ;
           func:score ?score .
    }
    FILTER(?score > 0.7)
}
```

Graph traversal stays in Rhea. Computation runs in the microservice. No data movement.

---

## Auto-generated docs and MCP

Google-style docstrings power both the README and an MCP resource:

```python
@mcp.resource("sparql://schema")
def sparql_schema() -> str:
    return ds.generate_docs()  # Markdown table: predicates, types, defaults

@mcp.prompt()
def sparql_assistant() -> str:
    return f"You are an expert in writing SPARQL queries...\n\n{ds.generate_docs()}"
```

LLMs can read the schema resource and write correct SPARQL queries against the service.

---

## Design Choices and Limitations

**What works well:**

- Pure Python: no RDF or SPARQL parsing for the implementer, type annotations and dataclasses drive all IRI and binding mapping, python defaults handle optional SPARQL inputs
- Works with federated `SERVICE` out of the box
- Enable building flexible Virtual Knowledge Graphs in a few lines of python

**Current limitations:**

- RDFLib uses a **global** custom evaluation registry: two `DatasetExt` instances in the same process share functions
- Custom function not supported for Oxigraph backend (possible)
- We force "conventions": python in snake_case, SPARQL in camelCase and PascalCase
- Performance: python-level evaluation per query, no query-planning awareness

---

## Thanks

**Install:**

```sh
uv add "rdflib-endpoint[web]"
pip install "rdflib-endpoint[web]"
```

**Library:** [github.com/vemonet/rdflib-endpoint](https://github.com/vemonet/rdflib-endpoint)

**Reference implementation:** [github.com/sib-swiss/mol-search-sparql-service](https://github.com/sib-swiss/mol-search-sparql-service)
