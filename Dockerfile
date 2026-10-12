# Build stage
FROM python:3.11-slim AS builder

WORKDIR /build
COPY . .
RUN pip install --no-cache-dir build && python -m build -w

# Final stage
FROM python:3.11-slim

RUN useradd -m -s /bin/bash bioseeker
USER bioseeker
WORKDIR /home/bioseeker

COPY --from=builder /build/dist/*.whl .
RUN pip install --no-cache-dir --user *.whl && rm *.whl

ENV PATH="/home/bioseeker/.local/bin:${PATH}"

ENTRYPOINT ["bioseeker"]
CMD ["--help"]
