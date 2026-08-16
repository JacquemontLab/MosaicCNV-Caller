FROM python:3.12-slim

ENV POLARS_VERSION=1.43.2

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*


RUN pip install --no-cache-dir "polars==${POLARS_VERSION}"