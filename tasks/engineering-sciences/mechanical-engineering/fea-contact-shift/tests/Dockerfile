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

RUN useradd --create-home --shell /usr/sbin/nologin runner

RUN mkdir -p /app
RUN mkdir -p /app/output
RUN mkdir -p /app/solver
RUN mkdir -p /app/data
RUN mkdir -p /logs/verifier

COPY . /tests/

# Practice packs are sealed into the verifier for instruction-alignment checks.
# /tests is root-only so the unprivileged analyzer child cannot read truth modules.
# test.sh / pytest still run as root in the verifier container.
RUN cp -a /tests/practice_data /app/data/practice \
 && chmod 755 /tests/test.sh \
 && chown -R root:root /tests /logs/verifier \
 && chmod -R u=rwX,go= /tests \
 && chmod 755 /tests/test.sh \
 && chmod 700 /logs/verifier \
 && chmod 755 /app /app/output /app/solver /app/data \
 && chmod -R a+rX /app/data/practice

WORKDIR /tests
