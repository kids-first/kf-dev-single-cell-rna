#!/usr/bin/env nextflow

include { MKREF } from './modules/local/mkref/main.nf'

workflow {
    main:
    // general
    fasta = params.fasta ? channel.fromPath(params.fasta) : error("Must provide a path to a FASTA file!")
    gtf = params.gtf ? channel.fromPath(params.gtf) : error("Must provide a path to a GTF file!")

    MKREF(
        fasta,
        gtf
    )
}