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
> The input TSV file is processed strictly by column order (headers are automatically detected and skipped). The engine will **stop and raise an error** if any row is structurally malformed (e.g. missing columns or invalid IRI syntax).
>
> 1.  **chem IRI**: IRI of the chemical. Must be explicitly wrapped in angle brackets (e.g., `<http://...>`, `<urn:...>`) or start with `http(s)://`.
> 2.  **SMILES**: SMILES string. If RDKit fails to parse the string into a valid molecule, the engine will **issue a warning and skip the row**.
> 3.  **db**: (Optional) Database source name. **Note on URIs vs Literals:** The engine preserves the exact string formatting provided. To specify a URI for accurate SPARQL matching later, you *must* wrap it in angle brackets (e.g., `<http://rhea>`). Otherwise, `Rhea` will be treated as a literal string.

Alternatively, you can dynamically fetch the compounds data on startup using a SPARQL query against a remote endpoint:

```bash
mol-search-sparql-service -s fetch_rhea.rq -e https://sparql.rhea-db.org/sparql
```

Other available optional flags include:
- `-t`, `--fingerprints`: Comma-separated list of fingerprint types to compute (e.g. `morgan_ecfp,pattern`). If omitted, all types are computed. The `pattern` fingerprint is **always computed** regardless of this option, as it is required for substructure search.
- `-p`, `--port`: Port to run the server on (default: `8010`).
- `-w`, `--workers`: Number of Uvicorn workers (default: `1`).
- `-d`, `--daemon`: Run the server in the background (logs to `server.log`).
- `-u`, `--public-url`: Public URL of the SPARQL endpoint exposed to clients (e.g. `https://api.example.com/sparql`). Used by the MCP server to advertise the correct endpoint address. Useful when running behind a reverse proxy. Defaults to `http://localhost:<port>/sparql`.

### 🌐 Default Endpoints

Once the service is running (default port 8010), you can access:

- **SPARQL Endpoint**: `http://localhost:8010/sparql` — Standard SPARQL Protocol endpoint for querying
- **MCP Server**: `http://localhost:8010/mcp` — Model Context Protocol server endpoint for LLM integration

### 📊 Query Examples

Once the service is running, you can query it using SPARQL. Here are some common examples:

**Similarity Search** — Find compounds similar to a query molecule:

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result ?score WHERE {
    [] a func:SimilaritySearch ;
        func:smiles "c1ccccc1" ;
        func:limit 5 ;
        func:result ?result ;
        func:score ?score .
}
```

**Substructure Search** — Find all compounds containing a specific substructure (e.g., benzene ring):

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result ?matchCount WHERE {
    [] a func:SubstructureSearch ;
        func:smart "c1ccccc1" ;
        func:limit 100 ;
        func:result ?result ;
        func:matchCount ?matchCount .
}
```

> [!NOTE]
> Substructure search requires the `pattern` fingerprint. It is **always computed automatically**, even when using the `-t/--fingerprints` option to restrict which fingerprints are loaded.

### 🧠 Memory Profile & Optimization

By default, the service computes **all 9 fingerprint types** (including chiral variants) for every loaded compound. Because they are hashed into `ExplicitBitVect` (2048-bits), each fingerprint consumes roughly **256 bytes per molecule**. This equates to a total footprint of approximately **~2.5 GB of RAM per 1,000,000 compounds**.

If you do not need all fingerprints, you can drastically reduce memory usage by explicitly passing the `--fingerprints` flag to specify only the ones you intend to search against.

> [!WARNING]
> **A Note on `atom_pair` and `topological_torsion`**: Previous versions of this service computed these fingerprints as exact, unhashed sparse vectors. For complex molecules, these sparse maps grew exponentially, causing memory usage to skyrocket. They are now forced into fixed-size 2048-bit dense vectors to cap memory and dramatically accelerate similarity speeds. While this is the industry standard for fast searching, minor bit-collisions may occur compared to exact sparse pairwise matching.


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


### `func:ListDatabases`

List available database names loaded in the service.

**IRI:** `urn:sparql-function:ListDatabases`

**Outputs:**

| Predicate | Type | Description |
|----------------------|------|-------------|
| `func:dbName` | `str` | The name or URI of the database stored in the service. |

**Example:**

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?dbName WHERE {
    [] a func:ListDatabases ;
        func:dbName ?dbName .
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
| `func:useChirality` | `bool` | `False` | If true, both tetrahedral (R/S) and double-bond (E/Z) stereochemistry are enforced during matching. Defaults to false. |

**Outputs:**

| Predicate | Type | Description |
|----------------------|------|-------------|
| `func:result` | `URIRef` | The URI of the matching compound. |
| `func:matchCount` | `int` | Number of matches found (1 if boolean match). |

**Example 1:**

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result ?matchCount WHERE {
    [] a func:SubstructureSearch ;
        func:smart "c1ccccc1" ;
        func:result ?result ;
        func:matchCount ?matchCount .
}
```

**Example 2:**

```sparql
PREFIX func: <urn:sparql-function:>
SELECT ?result ?matchCount WHERE {
    [] a func:SubstructureSearch ;
        func:smart "[C@@H](N)(O)F" ;
        func:useChirality true ;
        func:result ?result ;
        func:matchCount ?matchCount .
}
```

## 🫆 Available Fingerprint Types

| Key | Short Name | Description |
|-----|------------|-------------|
| `morgan_ecfp` | ECFP | Extended Connectivity Fingerprint (ECFP). Encodes atom-centered circular environments up to a given radius. Widely used for similarity search, clustering, and QSAR. |
| `morgan_ecfp_chiral` | ECFP_C | Extended Connectivity Fingerprint (ECFP). Encodes atom-centered circular environments up to a given radius. Widely used for similarity search, clustering, and QSAR. (Computed with stereochemistry enabled). |
| `morgan_fcfp` | FCFP | Functional-Class Fingerprint (FCFP). Morgan fingerprint using pharmacophoric atom features instead of exact atom types. |
| `morgan_fcfp_chiral` | FCFP_C | Functional-Class Fingerprint (FCFP). Morgan fingerprint using pharmacophoric atom features instead of exact atom types. (Computed with stereochemistry enabled). |
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

Generate slides:

```sh
cd slides
npm i
npm run dev
npm run build
```
