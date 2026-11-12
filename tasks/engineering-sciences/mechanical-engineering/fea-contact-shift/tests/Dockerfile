FROM python:3.12.7-slim-bookworm

ENV OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1

RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    scipy==1.14.1 \
    pytest==8.4.1 \
    pytest-json-ctrf==0.3.5

RUN useradd --create-home --shell /usr/sbin/nologin runner \
 && mkdir -p /app /app/output /app/solver /logs/verifier

COPY . /tests/

RUN chmod 755 /tests/test.sh \
 && chmod -R go-rwx /tests \
 && chmod 700 /logs/verifier \
 && chmod 755 /app /app/output /app/solver

WORKDIR /tests
