# ROADMAP: scRNA_seq_analysis

This document outlines the planned technical steps to develop a shareable and simple single-cell RNA-seq analysis pipeline, applicable to custom input data. This pipeline will facilitate reproducible and scalable scRNA-seq data processing and analysis.

## 🟡 Phase 1: Foundational Setup & Initial Data Ingestion

*   ✅ Task 1: Project Initialization & Version Control
    *   **Description**: Initialize the Git repository. Create a basic directory structure for Nextflow (`main.nf`, `nextflow.config`, `bin/` for scripts, `conf/` for profiles, `modules/` for Nextflow modules).
    *   **Validation**:
        *   Git repository is created (e.g., `git init`).
        *   Basic `main.nf` and `nextflow.config` files are created.
        *   Directory structure (`bin/`, `conf/`, `modules/`) is established.

*   ⚫ Task 2: Define Input Samplesheet Specification
    *   **Description**: Define the format for the input `samplesheet.csv` (e.g., columns: `sample_id`, `fastq_1`, `fastq_2`, `expected_cells`, `protocol_type`). Document this format within the README or a dedicated documentation file.
    *   **Validation**:
        *   A clear, documented specification for `samplesheet.csv` exists.
        *   An example `samplesheet.csv` file is created adhering to the specification.

*   ⚫ Task 3: Implement Basic FASTQ Input Reading in Nextflow
    *   **Description**: Develop a Nextflow process in `main.nf` or a dedicated module to:
        1.  Read the `samplesheet.csv`.
        2.  Validate the samplesheet structure and existence of FASTQ files.
        3.  Create Nextflow channels emitting tuples like `[meta_map, fastq_1, fastq_2]`.
    *   **Validation**:
        *   Nextflow pipeline successfully parses the example `samplesheet.csv` without errors.
        *   FASTQ file paths (and associated metadata) are correctly passed into Nextflow channels.
        *   The pipeline provides informative error messages for malformed samplesheets or missing files.

*   ⚫ Task 4: Initial Containerization Setup (Docker)
    *   **Description**: Configure Nextflow to use Docker. Create a simple Dockerfile for a basic bioinformatics tool that might be used early (e.g., FastQC for initial read quality assessment, or even just a Python environment for now).
    *   **Validation**:
        *   Nextflow successfully executes a simple process (e.g., listing files, running `echo`) using a Docker container defined by `docker.enabled = true` in `nextflow.config`.
        *   A basic Dockerfile (e.g., `environments/python/Dockerfile`) is created and builds successfully into an image.

*   ⚫ Task 5: Implement a Core Preprocessing/Quantification Module (e.g., Kallisto/BUStools)
    *   **Description**:
        1.  Choose one aligner/quantifier (e.g., Kallisto + BUStools as it's generally faster and lighter for initial development).
        2.  Create a Nextflow module (e.g., `modules/local/kallisto_bustools.nf`) to run the chosen tool.
        3.  This module should take FASTQ files and reference files (index, transcriptome, t2g map) as input.
        4.  It should output raw count matrices (e.g., in MTX format).
        5.  Containerize this module with Docker.
    *   **Validation**:
        *   The Nextflow module for quantification runs without errors on a small, well-defined test dataset (e.g., a subset of publicly available data).
        *   Expected output files (gene x cell matrix in MTX format, barcodes.tsv, genes.tsv) are generated in the correct structure.
        *   The process runs successfully within its Docker container.
        *   Parameters for reference files are configurable in `nextflow.config`.

## ⚫ Phase 2: Core Python-based Analysis Modules (Scanpy/AnnData)

*   ⚫ Task 1: Develop Python Script for Initial Data Loading & QC
    *   **Description**: Create a Python script (`bin/sc_qc.py`) using Scanpy/AnnData to:
        1.  Read the raw count matrix (MTX format) into an AnnData object.
        2.  Calculate standard QC metrics: n_genes_by_counts, total_counts, percentage_of_mitochondrial_reads (requires a gene list/GTF).
        3.  Allow configurable filtering based on these QC metrics.
        4.  Save the filtered AnnData object (`.h5ad` format).
        5.  Generate basic QC plots (e.g., violin plots of QC metrics before/after filtering).
    *   **Validation**:
        *   The Python script runs independently on example count data and produces a filtered `.h5ad` file.
        *   QC plots are generated as expected.
        *   QC thresholds are configurable via command-line arguments to the script.

*   ⚫ Task 2: Integrate QC Python Script into Nextflow
    *   **Description**:
        1.  Create a Nextflow module/process that takes the quantification output (MTX files) as input.
        2.  This process will execute the `sc_qc.py` script within its own container (ensure Python, Scanpy, Matplotlib are in the container).
        3.  Output the filtered `.h5ad` file and QC plots.
    *   **Validation**:
        *   The Nextflow process successfully runs the `sc_qc.py` script.
        *   The filtered `.h5ad` file and QC plots are correctly published to the results directory.
        *   The Python script's container is correctly defined and used.

*   ⚫ Task 3: Develop Python Script for Normalization & Feature Selection
    *   **Description**: Create a Python script (`bin/sc_normalize_select.py`) using Scanpy to:
        1.  Load the filtered `.h5ad` data.
        2.  Perform library size normalization (e.g., `sc.pp.normalize_total`).
        3.  Log-transform the data (e.g., `sc.pp.log1p`).
        4.  Identify highly variable genes (HVGs) (e.g., `sc.pp.highly_variable_genes`).
        5.  Subset the AnnData object to these HVGs.
        6.  Save the updated AnnData object.
    *   **Validation**:
        *   The Python script runs independently on a filtered `.h5ad` file.
        *   The output `.h5ad` file contains normalized data, log-transformed data, and HVG information.
        *   Parameters (e.g., for HVG selection) are configurable.

*   ⚫ Task 4: Integrate Normalization/Feature Selection Script into Nextflow
    *   **Description**: Create a Nextflow module/process to run the `sc_normalize_select.py` script, taking the QC'd `.h5ad` as input.
    *   **Validation**:
        *   The Nextflow process runs the script successfully.
        *   The updated `.h5ad` file (normalized, HVGs selected) is correctly published.

*   ⚫ Task 5: Develop Python Script for Dimensionality Reduction & Clustering
    *   **Description**: Create a Python script (`bin/sc_dimred_cluster.py`) using Scanpy to:
        1.  Load the normalized/feature-selected `.h5ad` data.
        2.  Perform Principal Component Analysis (PCA) (`sc.tl.pca`).
        3.  Compute a neighborhood graph (e.g., `sc.pp.neighbors` using PCA results).
        4.  Perform UMAP and/or t-SNE embedding (`sc.tl.umap`, `sc.tl.tsne`).
        5.  Perform cell clustering (e.g., Leiden algorithm: `sc.tl.leiden`).
        6.  Save the AnnData object with embeddings and cluster labels.
        7.  Generate UMAP/t-SNE plots colored by cluster labels and other relevant metadata.
    *   **Validation**:
        *   The Python script runs independently.
        *   The output `.h5ad` includes PCA, UMAP/t-SNE embeddings, and cluster labels.
        *   UMAP/t-SNE plots are generated.
        *   Key parameters (number of PCs, neighbors, resolution for clustering) are configurable.

*   ⚫ Task 6: Integrate Dimensionality Reduction/Clustering Script into Nextflow
    *   **Description**: Create a Nextflow module/process for this analysis step.
    *   **Validation**:
        *   The Nextflow process runs the script successfully.
        *   The final processed `.h5ad` and relevant plots are correctly published.

## ⚫ Phase 3: Advanced Analysis, HPC Support & Output Refinement

*   ⚫ Task 1: Develop Python Script for Marker Gene Identification
    *   **Description**: Create a Python script (`bin/sc_marker_genes.py`) using Scanpy to:
        1.  Load the clustered `.h5ad` data.
        2.  Find marker genes for each cluster (e.g., `sc.tl.rank_genes_groups` using t-test or Wilcoxon).
        3.  Save lists of marker genes (e.g., to CSV or TSV files).
        4.  Generate relevant summary plots (e.g., dot plots, heatmaps of top marker genes per cluster).
    *   **Validation**:
        *   The Python script runs independently on a clustered `.h5ad` file.
        *   Marker gene lists are generated for each cluster.
        *   Summary plots for marker genes are produced.

*   ⚫ Task 2: Integrate Marker Gene Script into Nextflow
    *   **Description**: Create a Nextflow module/process for marker gene identification.
    *   **Validation**:
        *   The Nextflow process runs the script successfully.
        *   Marker gene lists and plots are correctly published.

*   ⚫ Task 3: Implement Comprehensive Output Organization & Reporting
    *   **Description**:
        1.  Define a clear and user-friendly output directory structure using `publishDir` for all pipeline results (raw data links, quantification, QC reports, filtered data, processed data, plots, marker genes, etc.).
        2.  Consider implementing a MultiQC module to aggregate QC results from FastQC (if used) and custom Python script outputs.
    *   **Validation**:
        *   All pipeline outputs are organized logically and consistently in the specified output directory.
        *   A MultiQC report (if implemented) is generated and includes relevant metrics.
        *   The output structure is documented.

*   ⚫ Task 4: Apptainer/Singularity Container Support
    *   **Description**:
        1.  For each Docker container/environment defined, create an equivalent Apptainer/Singularity definition file.
        2.  Update `nextflow.config` to include a Singularity profile (`profiles.singularity`) that enables Singularity and specifies image locations if necessary.
    *   **Validation**:
        *   The entire pipeline runs successfully from start to finish using Singularity containers (`nextflow run . -profile singularity ...`).
        *   Singularity images are built/pulled correctly by Nextflow.

*   ⚫ Task 5: Refine Nextflow Configuration for HPC Execution
    *   **Description**:
        1.  Add or refine Nextflow profiles in `conf/` for common HPC schedulers (e.g., a generic `hpc` profile or specific ones like `slurm`).
        2.  Configure default resource requests (CPUs, memory, time) for different processes within these profiles.
        3.  Ensure users can easily override these defaults.
    *   **Validation**:
        *   The pipeline can be successfully submitted to and run on a test HPC environment (if accessible) using an HPC profile.
        *   Job resource requests are correctly interpreted by the HPC scheduler.

## ⚫ Phase 4: Testing, Documentation & Shareability

*   ⚫ Task 1: Create and Curate Small, Standardized Test Datasets
    *   **Description**: Prepare one or more small, publicly available (or mock) scRNA-seq datasets that can run quickly and test the full pipeline end-to-end, including different protocols if the pipeline aims to support them.
    *   **Validation**:
        *   Test datasets are available (e.g., in a `test_data/` directory or downloadable via a script).
        *   The pipeline successfully runs on these datasets from start to finish, producing expected outputs.
        *   The origin and characteristics of test data are documented.

*   ⚫ Task 2: Implement a Robust Pipeline Testing Strategy
    *   **Description**: Develop a strategy for testing pipeline integrity and reproducibility. This could involve:
        1.  Using `nf-test` framework for writing automated unit/integration tests for Nextflow modules and workflows.
        2.  Alternatively, simple shell scripts that run the pipeline on test data and check for the existence and basic integrity (e.g., non-zero file size, expected number of lines in a CSV) of key output files.
    *   **Validation**:
        *   Automated tests are implemented and can be run with a single command (e.g., `nf-test test` or `bash run_tests.sh`).
        *   Tests pass consistently on the CI (if set up) or locally for the main test datasets.

*   ⚫ Task 3: Write Comprehensive User Documentation
    *   **Description**: Consolidate and expand documentation into a user-friendly guide. This should include:
        1.  **README.md**: High-level overview, quick start, and links to detailed documentation.
        2.  **Installation**: How to install Nextflow, Docker/Singularity, and any other pipeline-specific dependencies.
        3.  **Input**: Detailed explanation of the `samplesheet.csv` and any other required inputs (e.g., reference genomes).
        4.  **Usage**: How to run the pipeline with examples for different profiles (local, Docker, Singularity, HPC).
        5.  **Parameters**: A comprehensive list of all configurable pipeline parameters, their defaults, and descriptions.
        6.  **Output**: Detailed explanation of the output directory structure and key files.
        7.  **Troubleshooting**: Common issues and how to resolve them.
    *   **Validation**:
        *   Documentation is clear, comprehensive, and accurate.
        *   A new user, unfamiliar with the pipeline, can successfully set up, run, and understand the results using only the provided documentation.

*   ⚫ Task 4: Finalize Codebase for Sharing & First Release
    *   **Description**:
        1.  Review all code (Nextflow scripts, Python scripts, container files) for clarity, consistency, and comments where necessary.
        2.  Ensure all dependencies are clearly specified (e.g., in container files, environment YAMLs if applicable).
        3.  Create a `LICENSE` file.
        4.  Create a `CHANGELOG.md` file.
        5.  Tag a version (e.g., `v1.0.0`) in Git.
        6.  If using GitHub/GitLab, create a release with release notes.
    *   **Validation**:
        *   The `README.md` is the definitive entry point and is fully up-to-date.
        *   The codebase is clean and well-organized.
        *   A `LICENSE` and `CHANGELOG.md` are present.
        *   The pipeline is tagged with a version number.
        *   The pipeline is ready for public sharing and use. 