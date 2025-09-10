process MULTI {
    label 'R8'
    container "migbro/cellranger:9.0.1"

    input:
        val(library_fastq_id)
        val(create_bam)
        path(reads)
        path(transcriptome)

    output:
    path("${library_fastq_id}/")

    script:
    def multi_csv = ""

    """
    cellranger telemetry disable;
    cellranger multi \\
    --id $library_fastq_id \\
    --csv $multi_csv 
    """
}