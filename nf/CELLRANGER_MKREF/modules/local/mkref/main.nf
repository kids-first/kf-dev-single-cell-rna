process MKREF {
    label 'CR'
    container "migbro/cellranger:9.0.1"

    input:
    path(fasta)
    path(gtf)

    output:
    path("*tar.gz")

    script:
    def out_genome_name = task.ext.prefix ?: "${gtf.baseName}"
    def gtf_args = task.ext.args
    """
    cellranger telemetry disable;
    cellranger mkgtf \\
    $gtf \\
    $gtf_args;
    echo "Filtered GTF created successfully"; \\
    cellranger mkref \\
    --genome=$out_genome_name \\
    --fasta=$fasta \\
    --genes=filtered.gtf; \\
    --disable-ui \\
    echo "Cell Ranger mkref finished successfully, packaging results";
    tar czvf ${out_genome_name}.tar.gz $out_genome_name/
    """
}