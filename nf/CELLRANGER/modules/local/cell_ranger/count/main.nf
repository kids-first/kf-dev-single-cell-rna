process COUNT {
    label 'CR'
    container "pgc-images.sbgenomics.com/d3b-bixu/cellranger:9.0.1"

    input:
        val(sample)
        path(reads)
        path(mates)
        path(transcriptome)
        path(index_files) // optional

    output:
    path("${sample}/outs/filtered_feature_bc_matrix.h5"), emit: filtered_matrix_h5
    path("${sample}/outs/raw_feature_bc_matrix.h5"), emit: raw_matrix_h5
    path("${sample}/outs/molecule_info.h5"), emit: molecule_info_h5
    path("${sample}/outs/analysis"), emit: analysis_dir, optional: true
    tuple path("${sample}/outs/possorted_genome_bam.bam"), path("${sample}/outs/possorted_genome_bam.bam.bai"), emit: bam_files, optional: true

    script:
    def link_reads = ""
    def flags = task.ext.args
    reads.eachWithIndex { read, index ->
        def mate = mates[index]
        link_reads += "cp -a $read fastqs/${sample}_S${index + 1}_L001_R1_001.fastq.gz && cp -a $mate fastqs/${sample}_S${index + 1}_L001_R2_001.fastq.gz;"
        if (index_files) {
            def index_file = index_files[index]
            link_reads += "cp -a $index_file fastqs/${sample}_S${index + 1}_L001_I1_001.fastq.gz;"
        }
    }
    """
    mkdir -p fastqs && $link_reads
    cellranger telemetry disable;
    cellranger count \\
    --disable-ui \\
    --id $sample \\
    --sample $sample \\
    --fastqs fastqs/ \\
    --transcriptome $transcriptome \\
    $flags;
    """
}