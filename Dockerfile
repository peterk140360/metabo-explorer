FROM python:3.10-slim

# System dependencies (needed esp. for rdkit)
RUN apt-get update && apt-get install -y \
    build-essential \
    ca-certificates \
    curl \
    libglib2.0-0 \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy dependency file first (better Docker caching)
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY src/ ./src/

# Expose Shiny port
EXPOSE 3838

# Environment defaults
ENV DATA_DIR=/data
ENV PYTHONUNBUFFERED=1

# Start the Shiny app
CMD ["shiny", "run", "--host", "0.0.0.0", "--port", "3838", "src/app.py"]
