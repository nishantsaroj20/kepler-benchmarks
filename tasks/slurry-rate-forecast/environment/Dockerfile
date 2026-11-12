FROM python:3.13-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        passwd \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 1000 --shell /usr/sbin/nologin runner

RUN pip install --no-cache-dir \
        iniconfig==2.3.1 \
        numpy==2.2.6 \
        packaging==26.3 \
        pluggy==1.7.0 \
        pygments==2.21.0 \
        pytest==8.4.1 \
        pytest-json-ctrf==0.3.5 \
        scipy==1.15.3

WORKDIR /app

COPY data_note.md /app/data_note.md
COPY data/ /app/data/

RUN mkdir -p /app/output /logs/verifier \
    && chmod 1777 /app /app/output /logs/verifier
