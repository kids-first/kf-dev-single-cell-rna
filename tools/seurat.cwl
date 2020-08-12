cwlVersion: v1.0
class: CommandLineTool
id: Seurat
doc: "Run Seurat analysis on 10x output"

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'satijalab/seurat:3.2.0'
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: seurat_analysis.R
        entry: |-
          #R script to run seurat analysis

          #this should be more manual, but will let us compare
          #the expected data to the actual

          #install packages if not installed.
          list.of.packages <- c("optparse")
          new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
          if(length(new.packages)) suppressMessages(install.packages(new.packages, "."))

          suppressMessages(library(optparse, lib.loc = "."))
          suppressMessages(library(patchwork))
          suppressMessages(library(dplyr))
          suppressMessages(library(Seurat))

          #command to check number of features(genes) and cells in seurat object
          #dim(GetAssayData(analysis))

          save_plot <- function(cmd, name) {
            #function to take in a command and a basename and output a plot file
            print(name)
            plot_file <- file.path(out_dir, paste(name, ".png", sep = ""))
            png(filename = plot_file, width = plot_size, height = plot_size)
            eval(parse(text = cmd))
          }

          #proccess inputs
          option_list <- list(
            make_option(
              opt_str = "--size",
              default = 500,
              type = "numeric",
              help = "Size for plots.",
            ),
            make_option(
              opt_str = "--name",
              default = "test",
              type = "character",
              help = "Project name."
            ),
            make_option(
              opt_str = "--out_size",
              default = 20,
              type = "numeric",
              help = "Number of genes to include in the cluster output file."
            ),
            make_option(
              opt_str = "--data",
              default = file.path(getwd(), "data"),
              type = "character",
              help = "Input data directory."
            ),
            make_option(
              opt_str = "--out",
              default = file.path(getwd(), "out"),
              type = "character",
              help = "Output directory path."
            )
          )

          #parse options
          opts <- parse_args(OptionParser(option_list = option_list))
          plot_size <- opts$size
          clust_size <- opts$out_size
          data_dir <- file.path(opts$data)
          out_dir <- file.path(opts$out)
          project_name <- opts$name

          #make output directory
          dir.create(out_dir, recursive = "true")

          ##load data
          analysis.data <- Read10X(data.dir = data_dir)
          analysis <- CreateSeuratObject(counts = analysis.data, project = project_name,
            min.cells = 3, min.features = 200)

          ##preprocess / qc
          #calculate % MT reads
          analysis[["percent.mt"]] <- PercentageFeatureSet(analysis, pattern = "^MT-")


          #Visualize QC metrics as a violin plot
          name <- "qc_violin"
          cmd <- 'VlnPlot(analysis, features = c("nFeature_RNA", "nCount_RNA",
            "percent.mt"), ncol = 3)'
          save_plot(cmd, name)

          #subset data with desired options
          analysis <- subset(analysis, subset = nFeature_RNA > 200 & nFeature_RNA < 2500
            & percent.mt < 20)

          #normalize data with selected type and scale factor
          analysis <- NormalizeData(analysis, normalization.method = "LogNormalize",
            scale.factor = 10000)

          #identify highly variable genes
          analysis <- FindVariableFeatures(analysis, selection.method = "vst",
            nfeatures = 2000)

          # Identify the 10 most highly variable genes
          top10 <- head(VariableFeatures(analysis), 10)

          #plot variable features with labels
          plot1 <- VariableFeaturePlot(analysis)
          plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
          name <- "features"
          cmd <- 'plot2'
          save_plot(cmd, name)

          ##pca

          #scale data
          all.genes <- rownames(analysis)
          analysis <- ScaleData(analysis, features = all.genes)

          #run PCA
          print("running PCA")
          analysis <- RunPCA(analysis, features = VariableFeatures(object = analysis))
          print("done with PCA")

          #output pca summary
          file <- file.path(out_dir, paste("pca_summary", ".txt", sep = ""))
          sink(file)
          print(analysis[["pca"]], dims = 1:7, nfeatures = clust_size)
          sink()

          #create a heat map for the first 10 PCs
          name <- "heat_map"
          cmd <- 'DimHeatmap(analysis, dims = 1:10, cells = 500, balanced = TRUE)'
          save_plot(cmd, name)

          #create an elbow plot for PCA
          name <- "elbow"
          cmd <- 'ElbowPlot(analysis)'
          save_plot(cmd, name)

          ##clustering
          analysis <- FindNeighbors(analysis, dims = 1:10)
          analysis <- FindClusters(analysis, resolution = 0.5)

          #run UMAP
          analysis <- RunUMAP(analysis, dims = 1:10)

          #plot UMAP
          name <- "umap"
          cmd <- 'DimPlot(analysis, reduction = "umap")'
          save_plot(cmd, name)

          #could also run tSNE at this point

          ##differential expression
          #find markers for every cluster compared to all remaining cells, report only the positive ones
          analysis.markers <- FindAllMarkers(analysis, only.pos = TRUE, min.pct = 0.25,
            logfc.threshold = 0.25)
          file <- file.path(out_dir, paste("cluster_markers", ".txt", sep = ""))
          cluster_markers <- analysis.markers %>% group_by(cluster) %>% top_n(n =
            clust_size, wt = avg_logFC)
          write.csv(cluster_markers, file = file)
          #print(cluster_markers)

          #the next few steps are just examples
          #eventually, we'll have to figure out what we would want

          #generate a heat map of the top 10 markers
          top10 <- analysis.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
          name <- "cluster_heat"
          cmd <- 'DoHeatmap(analysis, features = top10$gene) + NoLegend()'
          save_plot(cmd, name)

          ##save data object
          out_file <- file.path(out_dir, paste("analysis", ".rds", sep = ""))
          print("saving final data object")
          saveRDS(analysis, file = out_file)

        writable: false

baseCommand: [tar, -xaf]

arguments:
  - position: 1
    shellQuote: false
    valueFrom: >-
     $(inputs.data.path)
     && Rscript seurat_analysis.R --data $(inputs.data.nameroot.split('.')[0]) --out $(inputs.out)
     ${
       if (inputs.name != null){
         return "--name " + inputs.name;
       }
       else{
         return "";
       }
     }
     ${
       if (inputs.size != null){
         return "--size " + inputs.name;
       }
       else{
         return "";
       }
     }
     ${
       if (inputs.out_size != null){
         return "--out_size " + inputs.name;
       }
       else{
         return "";
       }
     }
     && tar -czf $(inputs.out).tar.gz $(inputs.out)

inputs:
  size: {type: int?, doc: "Plot size, plot area will be a square with side length of size"}
  name: {type: string?, doc: "Project string, used internall by Seurat"}
  out_size: {type: int?, doc: "Number of genes to include in the cluster output file."}
  data: {type: File, doc: "tarball of input data, folder must contain three files: barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz these files are the output of cellranger"}
  out: {type: string, doc: "name of output direcotry"}

outputs:
  out_file:
    type: File
    outputBinding:
      glob: $(inputs.out).tar.gz
    doc: "Tarball containing all files that were generated by Seurat"
