FROM --platform=linux/amd64 python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies (binary wheels only, no apt packages, no network at run time).
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir --only-binary=:all: -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12 -depth -type d -name __pycache__ -exec rm -rf {} + \
    && find /usr/local/lib/python3.12/site-packages -depth -type d \( -name tests -o -name test \) -exec rm -rf {} + \
    && rm -rf /usr/local/lib/python3.12/site-packages/pyarrow/include

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    HOME=/tmp

# Runs as uid 1000 on the platform; nothing here needs root. Nothing is written outside OUTPUT_DIR and /tmp.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.main"]
