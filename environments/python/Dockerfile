# Dockerfile for a basic Python environment for the scRNA-seq analysis pipeline

FROM python:3.10-slim

# Set a working directory
WORKDIR /opt/app

# Copy requirements file
# We copy it from the root of the project, assuming Nextflow will handle context
COPY requirements.txt .

# Install dependencies (currently requirements.txt is empty)
# RUN pip install --no-cache-dir -r requirements.txt

# Add a simple label
LABEL maintainer="scRNA_seq_analysis pipeline developer"
LABEL description="Basic Python environment for scRNA-seq Nextflow pipeline processes"

# The container will typically be run by Nextflow with a specific command,
# so no default CMD or ENTRYPOINT is strictly necessary here.
# For ad-hoc testing, one could add: CMD ["python"] 