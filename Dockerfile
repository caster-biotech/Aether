# Use official Python lightweight image
FROM python:3.11-slim

# Set environment variables to prevent Python from writing .pyc files and buffer outputs
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set the working directory in the container
WORKDIR /app

# Copy dependency specifications
COPY requirements.txt .

# Install production dependencies without caching to minimize image size
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source files and metadata
COPY src/ src/
COPY main.py .
COPY data/metadata/ data/metadata/

# Set the entrypoint to the main orchestration script
ENTRYPOINT ["python", "main.py"]
