FROM ghcr.io/astral-sh/uv:python3.13-trixie-slim
# https://docs.astral.sh/uv/guides/integration/docker

WORKDIR /app
COPY . /app/

RUN uv sync --frozen

ENV SERVER_PORT=8010
ENV WORKERS=1
ENV PYTHONUNBUFFERED='1'
# Path to SPARQL query file inside the container
ENV SPARQL_FILE='/app/fetch_rhea.rq'
# SPARQL endpoint to fetch compounds from
ENV SPARQL_ENDPOINT='https://sparql.rhea-db.org/sparql'

EXPOSE 8010
ENTRYPOINT ["sh", "-c", "uv run mol-search-sparql-service  -f compounds.tsv -p $SERVER_PORT -w $WORKERS"]
# ENTRYPOINT ["sh", "-c", "uv run mol-search-sparql-service -s $SPARQL_FILE -e $SPARQL_ENDPOINT -p $SERVER_PORT -w $WORKERS"]
