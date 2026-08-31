# syntax=docker/dockerfile:1

# -- Builder stage: install Python deps --
FROM python:3.11-slim AS builder

WORKDIR /app

# System deps for scipy/numpy compilation if needed (kept minimal)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir --prefix=/install -r requirements.txt

# -- Runtime stage --
FROM python:3.11-slim AS runtime

WORKDIR /app

# Create non-root user
RUN useradd --create-home --shell /bin/bash appuser

# Copy installed Python packages from builder
COPY --from=builder /install /usr/local

# Copy project code (only what the backend needs at runtime)
COPY backend/ ./backend/
COPY data_pipeline/ ./data_pipeline/
COPY model/ ./model/
COPY requirements.txt .

# Ensure the app package is importable (WORKDIR is on PYTHONPATH, but be explicit)
ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1

# Release gate env vars are NOT set here — fail-closed by default.
# Override at runtime via docker-compose or `docker run -e`:
#   MODEL_RELEASE_APPROVED=true
#   APPROVED_ARTIFACT_REVISION=<revision>
#   MODEL_ARTIFACT_PATH=/app/data/graph.json
# The backend's model_store.py enforces fail-closed logic; this image does not weaken it.

EXPOSE 8000

USER appuser

# Healthcheck hits the liveness endpoint (always 200 regardless of release gate)
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health').read()" || exit 1

CMD ["uvicorn", "backend.app:app", "--host", "0.0.0.0", "--port", "8000"]
