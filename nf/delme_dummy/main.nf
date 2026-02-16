process list_stuff {
    label 'process_single'
    input:
    path(inputs)
    output:
    stdout
    script:
    """
    ls -lH ${inputs.join(' ')}
    """
}
workflow {
    main:

    manifest = file(params.manifest)
    inputs = channel.fromPath(manifest)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            [row.file_name, row.dir_name]
                .findAll { it && it.trim() != '' }
                .collect { file(it) }
        }
        .flatten()
    list_stuff(inputs.collect())
}