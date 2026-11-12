FROM python:3.12.7-slim-bookworm

ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1

RUN pip install --no-cache-dir numpy==1.26.4 scipy==1.14.1

RUN mkdir -p /app/solver /app/output /app/data

COPY data/ /app/data/

WORKDIR /app
