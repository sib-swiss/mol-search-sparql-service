# 🫆 Chemistry Search SPARQL Service

A service for finding similar chemicals and substructures using molecular fingerprints. It provides both a SPARQL endpoint and an MCP (Model Context Protocol) server interface.

Built using [**RDKit**](https://www.rdkit.org/) and [**rdflib-endpoint**](https://github.com/vemonet/rdflib-endpoint).

## 🧩 Features

- **SPARQL endpoint**: Query for similarity and substructure matches using standard SPARQL syntax. You can also query this service from another SPARQL endpoint using the `SERVICE` clause:

  ```sparql
  PREFIX func: <urn:sparql-function:>
  SELECT ?result ?score WHERE {
    SERVICE <http://localhost:8000/sparql> {
      [] a func:SimilaritySearch ;
         func:smiles "c1ccccc1" ;
         func:result ?result ;
         func:score ?score .
    }
  }
  ```

- **MCP Server**: Exposes SPARQL functions documentation through resource and prompt to help LLMs to perform searches (experimental).

## 📥 Installation

You can install this tool directly from GitHub using `uv`:

```bash
uv tool install git+https://github.com/sib-swiss/mol-search-sparql-service.git
```

## 🚀 Usage

Start the service by pointing it to a compounds data file (TSV):

```bash
mol-search-sparql-service -f compounds.tsv
```

> [!NOTE]
>
> The input TSV file should contain at least:
>
> -   `?chem`: IRI of the chemical
> -   `?smiles`: SMILES string
> -   `?db`: optional database source URI

Alternatively, you can dynamically fetch the compounds data on startup using a SPARQL query against a remote endpoint:

```bash
mol-search-sparql-service -s fetch_rhea.rq -e https://sparql.rhea-db.org/sparql
```

Other available optional flags include `-w` (`--workers`) to deploy multiple Uvicorn worker processes (default 1) and `-d` (`--daemon`) to run the server in the background and write stdout/stderr to `server.log`. Port defaults to `8010` if `-p` is omitted.

<!-- AUTOGEN_DOCS_START -->

## 📖 Functions

### `func:ListFingerprints`

List available fingerprint types.

**IRI:** `urn:sparql-function:ListFingerprints`

**Outputs:**

| Predicate | Type | Description |
|----------------------|------|-------------|
| `func:fpType` | `str` | Identifier key for the fingerprint type (e.g. `morgan_ecfp`). |
| `func:description` | `str` | Human readable description of the fingerprint. |
| `func:mechanism` | `str` | Explainability mechanism / how bits map to substructures. |
| `func:shortName` | `str` | Short display name for the fingerprint (e.g. ECFP, MACCS). |

**Example:**

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?fpType ?description ?shortName WHERE {
    [] a func:ListFingerprints ;
        func:fpType ?fpType ;
        func:description ?description ;
        func:shortName ?shortName .
}
```


### `func:SimilaritySearch`

Perform similarity search using precomputed fingerprints.

**IRI:** `urn:sparql-function:SimilaritySearch`

**Inputs:**

| Predicate | Type | Default | Description |
|-----------------|------|---------|-------------|
| `func:smiles` | `str` | *required* | Query SMILES string. |
| `func:limit` | `int` | `10` | Maximum number of results to return. |
| `func:dbNames` | `UnionType[str, NoneType]` | `None` | Optional database name to filter results. |
| `func:fpType` | `str` | `'morgan_ecfp'` | Fingerprint type key to use. |
| `func:useChirality` | `bool` | `False` | Whether to respect chirality when computing fingerprints. |
| `func:minScore` | `float` | `0.0` | Minimum similarity score threshold (0.0 - 1.0). |

**Outputs:**

| Predicate | Type | Description |
|----------------------|------|-------------|
| `func:result` | `URIRef` | The URI of the matching compound. |
| `func:score` | `float` | Tanimoto similarity score (0-1). |

**Example:**

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


### `func:SubstructureSearch`

Perform substructure search using a SMARTS/SMILES pattern.

**IRI:** `urn:sparql-function:SubstructureSearch`

**Inputs:**

| Predicate | Type | Default | Description |
|-----------------|------|---------|-------------|
| `func:smart` | `str` | *required* | Query SMARTS or SMILES pattern to match. |
| `func:limit` | `int` | `100` | Maximum number of results to return (default: 100). |
| `func:dbNames` | `UnionType[str, NoneType]` | `None` | Optional database name to limit the search. |
| `func:useChirality` | `bool` | `False` | Whether to respect chirality for matching. |
| `func:minMatchCount` | `int` | `1` | Minimum number of substructure matches required. |

**Outputs:**

| Predicate | Type | Description |
|----------------------|------|-------------|
| `func:result` | `URIRef` | The URI of the matching compound. |
| `func:matchCount` | `int` | Number of matches found (1 if boolean match). |

**Example:**

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result ?matchCount WHERE {
    [] a func:SubstructureSearch ;
        func:smart "c1ccccc1" ;
        func:result ?result ;
        func:matchCount ?matchCount .
}
```

## 🫆 Available Fingerprint Types

| Key | Short Name | Description |
|-----|------------|-------------|
| `morgan_ecfp` | ECFP | Extended Connectivity Fingerprint (ECFP). Encodes atom-centered circular environments up to a given radius. Widely used for similarity search, clustering, and QSAR. |
| `morgan_fcfp` | FCFP | Functional-Class Fingerprint (FCFP). Morgan fingerprint using pharmacophoric atom features instead of exact atom types. |
| `rdk_topological` | RDK | RDKit topological (path-based) fingerprint. Encodes linear bond paths similar to Daylight fingerprints. |
| `atom_pair` | AP | Atom Pair fingerprint. Encodes pairs of atoms along with their topological distance. |
| `topological_torsion` | TT | Topological Torsion fingerprint. Encodes sequences of four bonded atoms. |
| `maccs` | MACCS | MACCS structural keys (166 bits). Each bit corresponds to a predefined chemical pattern. |
| `pattern` | Pattern | RDKit Pattern fingerprint. Designed for substructure screening. |

<!-- AUTOGEN_DOCS_END -->

## 🛠️ Development

Run tests:

```sh
uv run pytest
```

Format and lint:

```sh
uvx ruff format && uvx ruff check --fix
```

Auto-generate docs from functions docstrings and update the `README.md`:

```sh
uv run src/mol_search_sparql_service/gen_docs.py
```

To release a new version run the release script providing the version bump: `fix`, `minor`, or `major`:

```sh
.github/release.sh fix
```

