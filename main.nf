#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters
// params.input is now expected to be defined in nextflow.config or command line
// params.outdir is also expected from nextflow.config or command line

log.info """
R U N N I N G - P I P E L I N E
===================================
Input Samplesheet: ${params.input}
Output Directory:  ${params.outdir}
===================================
"""

// Channel for reading the samplesheet
Channel
    .fromPath(params.input)
    .ifEmpty { error "Cannot find samplesheet file: ${params.input}" }
    .splitCsv(header:true, sep:',')
    .map { row ->
        // Validate mandatory columns
        if (!row.sample || !row.fastq_1 || !row.protocol) {
            error "Samplesheet CSV '${params.input}' is missing one or more mandatory columns (sample, fastq_1, protocol) in row: ${row}"
        }

        // Create meta map
        def meta = [:]
        meta.id = row.sample.trim()
        meta.protocol = row.protocol.trim()
        if (row.expected_cells && !row.expected_cells.trim().isEmpty()) {
            try {
                meta.expected_cells = row.expected_cells.trim().toInteger()
            } catch (NumberFormatException e) {
                log.warn "Warning: Could not parse 'expected_cells' for sample '${meta.id}' (value: '${row.expected_cells}') as an integer. It will be ignored for this sample."
                meta.expected_cells = null
            }
        } else {
            meta.expected_cells = null
        }

        def fq1_path = row.fastq_1.trim()
        if (fq1_path.isEmpty()) {
            error "FASTQ_1 path is empty for sample '${meta.id}' in samplesheet '${params.input}'."
        }
        def fq1_file = file(fq1_path)
        if (!fq1_file.exists()) {
            error "FASTQ file for read 1 not found for sample '${meta.id}': ${fq1_path} (resolved to: ${fq1_file.getAbsolutePath()})"
        }

        def fq2_file = null
        if (row.fastq_2 && !row.fastq_2.trim().isEmpty()) {
            def fq2_path = row.fastq_2.trim()
            fq2_file = file(fq2_path)
            if (!fq2_file.exists()) {
                error "FASTQ file for read 2 not found for sample '${meta.id}': ${fq2_path} (resolved to: ${fq2_file.getAbsolutePath()})"
            }
        }
        def fq2_to_emit = fq2_file ?: [] // Use empty list if fq2_file is null
        return [ meta, fq1_file, fq2_to_emit ]
    }
    .set { reads_ch }

// Process to demonstrate Docker execution
process SIMPLE_DOCKER_ECHO {
    tag "echo_${meta.id}"
    publishDir "${params.outdir}/simple_docker_echo", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(fq1), path(fq2)

    output:
    path "${meta.id}_echo.txt"

    container 'docker://alpine:latest' // Explicitly use docker URI for Apptainer

    script:
    """
    echo "Hello from Docker container for sample ${meta.id}!"
    echo "FASTQ1: ${fq1.name}"
    if [ -n "${fq2}" ]; then // Check if fq2 is a non-empty string (i.e., a file was provided)
        echo "FASTQ2: ${fq2.name}"
    else
        echo "FASTQ2: N/A"
    fi
    echo "Outputting to ${meta.id}_echo.txt"
    echo "Hello from Docker container for sample ${meta.id}! FASTQ1: ${fq1.name}" > ${meta.id}_echo.txt
    """
}

// Workflow definition
workflow {
    main:
        // Original view for samplesheet parsing (can be commented out later)
        reads_ch.view { meta, fq1, fq2 ->
            def fq2_display = (fq2 instanceof List && fq2.isEmpty()) || fq2 == null ? 'N/A' : fq2.name
            def expected_cells_display = meta.expected_cells != null ? meta.expected_cells : 'N/A'
            "Sample (view): ${meta.id}, Protocol: ${meta.protocol}, Expected Cells: ${expected_cells_display}, FASTQ1: ${fq1.name}, FASTQ2: ${fq2_display}"
        }

        // Call the Dockerized echo process
        SIMPLE_DOCKER_ECHO(reads_ch)

        // View the output of the Dockerized process
        SIMPLE_DOCKER_ECHO.out.view { file_path ->
            "Docker echo output file: ${file_path}"
        }
}

workflow.onComplete {
    log.info ( workflow.success ? """
        Pipeline complete!
        Output directory: ${params.outdir}
        """ : """
        Pipeline FAILED!
        See log file for details (e.g., .nextflow.log or check SLURM/LSF logs if applicable)
        """
    )
} 