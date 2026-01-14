# Base image
FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY . .

# Run build script
RUN ./build.sh

# Set entrypoint
ENTRYPOINT ["bioseeker"]
CMD ["--help"]
