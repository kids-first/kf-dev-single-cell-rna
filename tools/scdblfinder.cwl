cwlVersion: v1.2
class: CommandLineTool
id: scdblfinder-rmd
doc: "Run custom QC on 10X output"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/d3b-bixu/scdblfinder:1.12.0"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: scDblFinder.Rmd
        entry:
          $include: ../scripts/scDblFinder.Rmd
  - class: ResourceRequirement
    coresMin: $(inputs.cpus)
    ramMin: ${ return inputs.ram * 1000 }

baseCommand: [Rscript, -e]
arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
      "rmarkdown::render('scDblFinder.Rmd', clean = TRUE,
            params=list(results_dir='.', 
                        data_path='$(inputs.seurat_raw_object.path)', 
                        sample_name='$(inputs.sample_name)',
                        cpus=$(inputs.cpus),
                        ram=$(inputs.ram)
                        ))"
inputs:
  seurat_raw_object: { type: File, doc: "Seurat RDS from align/qc workflow" }
  sample_name: { type: string }
  cpus: { type: 'int?', doc: "Num CPUs for Seurat parallization",
    default: 8 }
  ram: { type: 'int?', doc: "Ram in GB for task",
    default: 16 }

outputs:
  result_dir:
    type: Directory
    outputBinding:
      loadListing: deep_listing
      glob: 'scDblFinder-$(inputs.sample_name)'
    doc: "Dir containing all files that were generated by the notebook"
