FROM ghcr.io/astral-sh/uv:python3.13-trixie-slim
# https://docs.astral.sh/uv/guides/integration/docker

WORKDIR /app

# RDKit's drawing module (rdMolDraw2D) needs X11 shared libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 libxext6 libsm6 libexpat1 \
    && rm -rf /var/lib/apt/lists/*

COPY . /app/

RUN uv sync --frozen

ENV SERVER_PORT=8010
ENV WORKERS=1
ENV PYTHONUNBUFFERED='1'
# Path to SPARQL query file inside the container
ENV SPARQL_FILE='/app/fetch_rhea.rq'
# SPARQL endpoint to fetch compounds from
ENV SPARQL_ENDPOINT='https://sparql.rhea-db.org/sparql'

# Precompute fingerprints from the bundled compounds.tsv so the server starts
# fast (no recompute on cold start). The cache is keyed on file content + fp
# types, so the runtime entrypoint below reuses it directly.
RUN uv run python -c "from mol_search_sparql_service.rdkit_fingerprints import engine; engine.load_file('compounds.tsv', cache_dir='.mol-search-service')"

# Run as a non-root user (Hugging Face Spaces convention). Give it ownership of
# /app so the precomputed cache and venv remain readable/writable at runtime.
RUN useradd -m -u 1000 user && chown -R user:user /app
USER user
ENV HOME=/home/user

EXPOSE 8010
ENTRYPOINT ["sh", "-c", "uv run mol-search-sparql-service -f compounds.tsv -p $SERVER_PORT -w $WORKERS"]
# ENTRYPOINT ["sh", "-c", "uv run mol-search-sparql-service  -s fetch_rhea.rq -e https://sparql.rhea-db.org/sparql -p $SERVER_PORT -w $WORKERS"]
# ENTRYPOINT ["sh", "-c", "uv run mol-search-sparql-service -s $SPARQL_FILE -e $SPARQL_ENDPOINT -p $SERVER_PORT -w $WORKERS"]
