# 🫆 Chemistry Search SPARQL Service

A service for finding similar chemicals and substructures using molecular fingerprints. It provides both a SPARQL endpoint and an MCP (Model Context Protocol) server interface.

Built using [**RDKit**](https://www.rdkit.org/) and [**rdflib-endpoint**](https://github.com/vemonet/rdflib-endpoint).

## 🧩 Features

- **SPARQL endpoint**: Query for similarity and substructure matches using standard SPARQL syntax. You can also query this service from another SPARQL endpoint using the `SERVICE` clause:

  ```sparql
  PREFIX func: <urn:sparql-function:>
  SELECT ?result ?score WHERE {
    SERVICE <http://localhost:8010/sparql> {
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

> [!IMPORTANT]
>
> The input TSV file is processed strictly by column order (headers are automatically detected and skipped). The engine will **stop and raise an error** if any row is malformed.
>
> 1.  **chem IRI**: IRI of the chemical. Must be explicitly wrapped in angle brackets (e.g., `<http://...>`, `<urn:...>`) or start with `http(s)://`.
> 2.  **SMILES**: SMILES string. Must be successfully parseable into a valid molecule by RDKit.
> 3.  **db**: (Optional) Database source name. **Note on URIs vs Literals:** The engine preserves the exact string formatting provided. To specify a URI for accurate SPARQL matching later, you *must* wrap it in angle brackets (e.g., `<http://rhea>`). Otherwise, `Rhea` will be treated as a literal string.

Alternatively, you can dynamically fetch the compounds data on startup using a SPARQL query against a remote endpoint:

```bash
mol-search-sparql-service -s fetch_rhea.rq -e https://sparql.rhea-db.org/sparql
```

Other available optional flags include:
- `-p`, `--port`: Port to run the server on (default: `8010`).
- `-w`, `--workers`: Number of Uvicorn workers (default: `1`).
- `-d`, `--daemon`: Run the server in the background (logs to `server.log`).


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
| `morgan_ecfp_chiral` | ECFP_C | Extended Connectivity Fingerprint (ECFP). Encodes atom-centered circular environments up to a given radius. Widely used for similarity search, clustering, and QSAR. (Computed with stereochemistry enabled). |
| `morgan_fcfp_chiral` | FCFP_C | Functional-Class Fingerprint (FCFP). Morgan fingerprint using pharmacophoric atom features instead of exact atom types. (Computed with stereochemistry enabled). |

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

Generate slides:

```sh
cd slides
npm i
npm run dev
npm run build
```
