#!/usr/bin/env nextflow

include { MKREF } from './modules/local/mkref/main.nf'

workflow {
    main:
    // general
    fasta = params.fasta ? Channel.fromPath(params.fasta) : error("Must provide a path to a FASTA file!")
    gtf = params.gtf ? Channel.fromPath(params.gtf) : error("Must provide a path to a GTF file!")

    MKREF(
        fasta,
        gtf
    )
}