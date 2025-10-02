process MKREF {
    label 'CR'
    container "migbro/cellranger:9.0.1"

    input:
    path(fasta)
    path(gtf)
    val(out_genome_name)
    val(attribute_filter)

    output:
    path("${out_genome_name}.tar.gz")

    script:
    """
    cellranger telemetry disable;
    cellranger mkgtf \\
    $gtf \\
    filtered.gtf \\
    --attribute=$attribute_filter;
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