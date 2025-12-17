process COUNT {
    label 'CR'
    container "pgc-images.sbgenomics.com/d3b-bixu/cellranger:9.0.1"

    input:
        val(sample)
        val(create_bam)
        path(reads)
        path(mates)
        path(transcriptome)
        path(index_files) // optional

    output:
    path("${sample}.tar.gz")

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
    $flags; \\
    echo "Cell Ranger count finished successfully, packaging results";
    tar czvf ${sample}.tar.gz $sample/outs
    """
}