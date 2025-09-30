#!/usr/bin/env nextflow

include { MKREF } from '../CELLRANGER/modules/local/cell_ranger/mkref/main.nf'

workflow {
    main:
    // general
    fasta = params.fasta ? file(params.fasta) : error("Must provide a path to a FASTA file!")
    gtf = params.gtf ? file(params.gtf) : error("Must provide a path to a GTF file!")
    out_genome_name = params.out_genome_name ? params.out_genome_name : error("Must provide a name for the output genome!")
    attribute_filter = params.attribute_filter ? params.attribute_filter : "gene_biotype:protein_coding"

    MKREF(
        fasta,
        gtf,
        out_genome_name,
        attribute_filter
    )
}