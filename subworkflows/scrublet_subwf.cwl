cwlVersion: v1.2
class: Workflow
id: scrublet_subwf

doc: "Subworkflow to run scrublet on matrix files in 10x output directory."

requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  count_dir: {type: Directory, doc: "Directory containing raw and filtered mtx files produced by Cell Ranger count"}
  output_basename: {type: string, doc: "Output files basename"}
  expected_doublet_rate: {type: 'float?'}
  doublet_score_threshold: {type: 'float?'}
  count_min: {type: 'int?'}
  cell_min: {type: 'int?'}
  min_gene_variability_pctl: {type: 'int?'}
  n_prin_comps: {type: 'int?'}
  ram: {type: 'int?', default: 16, doc: "In GB"}
  cpus: {type: 'int?', default: 1, doc: "Number of CPUs to request"}

outputs:
  out_dir: {type: Directory, outputSource: add_filtered_matrices/out_dir}
  filtered_doublet_histogram: {type: File, outputSource: run_scrublet_filtered/score_histogram}
  raw_doublet_histogram: {type: File, outputSource: run_scrublet_raw/score_histogram}

steps:

  find_matrix_files:
    run: ../tools/find_matrix_files.cwl
    in:
      count_dir: count_dir
    out: [filtered_matrix, raw_matrix]

  run_scrublet_filtered:
    run: ../tools/scrublet.cwl
    in:
      input_matrix: find_matrix_files/filtered_matrix
      output_basename: output_basename
      expected_doublet_rate: expected_doublet_rate
      doublet_score_threshold: doublet_score_threshold
      count_min: count_min
      cell_min: cell_min
      min_gene_variability_pctl: min_gene_variability_pctl
      n_prin_comps: n_prin_comps
      ram: ram
      cpus: cpus
    out: [score_histogram, doublet_removed_matrix]

  run_scrublet_raw:
    run: ../tools/scrublet.cwl
    in:
      input_matrix: find_matrix_files/filtered_matrix
      output_basename: output_basename
      expected_doublet_rate: expected_doublet_rate
      doublet_score_threshold: doublet_score_threshold
      count_min: count_min
      cell_min: cell_min
      min_gene_variability_pctl: min_gene_variability_pctl
      n_prin_comps: n_prin_comps
      ram: ram
      cpus: cpus
    out: [score_histogram, doublet_removed_matrix]

  add_filtered_matrices:
    run: ../tools/add_filtered_matrices.cwl
    in:
      count_dir: count_dir
      filtered_matrix: run_scrublet_filtered/doublet_removed_matrix
      raw_matrix: run_scrublet_raw/doublet_removed_matrix
    out: [out_dir]
