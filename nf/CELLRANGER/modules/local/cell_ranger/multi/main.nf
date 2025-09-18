process MULTI {
    label 'R8'
    container "migbro/cellranger:9.0.1"

    input:
        val(library_fastq_id)
        val(create_bam)
        path(reads, name: "FASTQS/")
        path(mates, name: "FASTQS/")
        path(transcriptome)
        val(feature_types)
        path(sample_sheet)
        path(probe_set_csv)

    output:
    path("${library_fastq_id}/outs/multi/"), emit: multi_out
    path("${library_fastq_id}/outs/per_sample_outs"), emit: per_sample
    path("*multi_config.csv"), emit: config

    script:
    def multi_config = "${library_fastq_id}.multi_config.csv"

    """
    cellranger telemetry disable;
    echo -e "[gene-expression]
    reference,\$PWD/$transcriptome
    probe-set,\$PWD/$probe_set_csv
    create-bam,$create_bam

    [libraries]
    fastq_id,fastqs,feature_types
    $library_fastq_id,\$PWD/FASTQS/,$feature_types
    
    [samples]" > $multi_config && \\
    cat $sample_sheet >> $multi_config && \\
    cellranger multi \\
    --id $library_fastq_id \\
    --csv $multi_config 
    """
}