FROM python:3.13-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        iniconfig==2.3.1 \
        packaging==26.3 \
        pluggy==1.7.0 \
        pygments==2.21.0 \
        pytest==8.4.1 \
        pytest-json-ctrf==0.3.5

COPY . /tests/

RUN mkdir -p /app /app/output /logs/verifier
