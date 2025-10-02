#!/usr/bin/env nextflow

include { MKREF } from './modules/local/mkref/main.nf'

workflow {
    main:
    // general
    fasta = params.fasta ? Channel.fromPath(params.fasta) : error("Must provide a path to a FASTA file!")
    gtf = params.gtf ? Channel.fromPath(params.gtf) : error("Must provide a path to a GTF file!")
    out_genome_name = params.out_genome_name ? Channel.value(params.out_genome_name) : error("Must provide a name for the output genome!")
    attribute_filter = params.attribute_filter ? Channel.value(params.attribute_filter) : error("Must provide an attribute_filter")

    MKREF(
        fasta,
        gtf,
        out_genome_name,
        attribute_filter
    )
}