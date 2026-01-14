process MULTI {
    label 'CR'
    container "pgc-images.sbgenomics.com/d3b-bixu/cellranger:9.0.1"

    input:
        val(library_fastq_id)
        path(reads, name: "FASTQS/")
        path(mates, name: "FASTQS/")
        path(transcriptome)
        val(feature_types)
        path(sample_sheet)
        path(probe_set_csv)

    output:
    path("${library_fastq_id}/outs/per_sample_outs/*/count/analysis"), emit: multi_analysis, optional: true
    path("${library_fastq_id}/outs/per_sample_outs/*/count/*.h5"), emit: multi_h5
    path("${library_fastq_id}/outs/per_sample_outs/*/metrics_summary.csv"), emit: multi_metrics
    tuple path("${library_fastq_id}/outs/per_sample_outs/*/count/*.bam"), path("${library_fastq_id}/outs/per_sample_outs/*/count/*.bai"), emit: multi_bam, optional: true
    path("*multi_config.csv"), emit: config

    script:
    def multi_config = "${library_fastq_id}.multi_config.csv"
    def flags = task.ext.args

    """
    cellranger telemetry disable;
    echo -e "[gene-expression]
    reference,\$PWD/$transcriptome
    probe-set,\$PWD/$probe_set_csv
    $flags

    [libraries]
    fastq_id,fastqs,feature_types
    $library_fastq_id,\$PWD/FASTQS/,$feature_types
    
    [samples]" > $multi_config; \\
    cat $sample_sheet >> $multi_config; \\
    cat $multi_config; \\
    cellranger multi \\
    --disable-ui \\
    --id $library_fastq_id \\
    --csv $multi_config
    """
}