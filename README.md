# Chemistry Search SPARQL Service

A SPARQL service for finding similar chemicals and substructures using molecular fingerprints. Built using **RDKit**.

## Installation

You can install this tool directly from GitHub using `uv`:

```bash
uv tool install git+https://github.com/sib-swiss/mol-search-sparql-service.git
```

## Usage

Start the service by pointing it to a compounds data file (TSV):

```bash
mol-search-sparql-service -f compounds.tsv -p 8000
```

Alternatively, you can dynamically fetch the compounds data on startup using a SPARQL query against a remote endpoint:

```bash
mol-search-sparql-service -s fetch_rhea.rq -e https://sparql.rhea-db.org/sparql -p 8000
```

Other available optional flags include `-w` (`--workers`) to deploy multiple Uvicorn worker processes (default 1) and `-d` (`--daemon`) to run the server in the background and write stdout/stderr to `server.log`. Port defaults to `8010` if `-p` is omitted.

## Features

-   **SPARQL Endpoint**: Query for similarity and substructure matches using standard SPARQL syntax.
-   **MCP Server**: Exposes tools for LLMs to perform searches.
-   **Introspection**: Discover available fingerprint types and schema details.

### Federated Query (SERVICE clause)

You can also query this service from another SPARQL endpoint using the `SERVICE` clause:

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

### Data Format
The input TSV file should contain at least:
-   `?chem`: URI of the chemical
-   `?smiles`: SMILES string
-   `?db`: (Optional) Database source URI

## Development

Run tests:

```sh
uv run pytest
```

