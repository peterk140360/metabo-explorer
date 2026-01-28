# =========================
# Base image
# =========================
FROM python:3.12-slim-bookworm

# =========================
# System dependencies
# =========================
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    libglib2.0-0 \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

# =========================
# App directory
# =========================
WORKDIR /app/src

# =========================
# Python dependencies
# =========================
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r /app/requirements.txt

# =========================
# Copy application code
# =========================
COPY src/ /app/src/

# =========================
# Runtime configuration
# =========================
ENV DATA_DIR=/app/data
ENV PYTHONUNBUFFERED=1

# =========================
# Shiny port
# =========================
EXPOSE 3838

# =========================
# Start Shiny app
# =========================
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "3838", "app.py"]
