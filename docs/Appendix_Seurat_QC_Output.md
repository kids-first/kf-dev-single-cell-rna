# D3b Seurat QC Output File
Created by @AntoniaChroni, this object is the result of an R Markdown notebook, which in summary does the following:
 - Summarize Feature (cell) count, RNA (gene/cell) count, and percent mitochondria via violin plot
 - Apply a basics min-max gene count filter as well as max percent mitchondria
 - Plot most expressed genes
 - Plot post-filter violin plots
 - Normalize and plot PCA
 - Plot various UMAPs of features by attribute

Below is the `str()` output of the resulting Seurat Object rds file:

```
Formal class 'Seurat' [package "SeuratObject"] with 13 slots
  ..@ assays      :List of 1
  .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  .. .. .. .. .. ..@ i       : int [1:5954964] 52 58 71 80 88 103 107 176 197 201 ...
  .. .. .. .. .. ..@ p       : int [1:3224] 0 1544 3756 5145 6744 9504 10322 11496 13834 17232 ...
  .. .. .. .. .. ..@ Dim     : int [1:2] 33425 3223
  .. .. .. .. .. ..@ Dimnames:List of 2
  .. .. .. .. .. .. ..$ : chr [1:33425] "DDX11L1" "WASH7P" "MIR1302-2HG" "ENSG00000238009" ...
  .. .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..@ x       : num [1:5954964] 1 3 1 5 1 4 5 13 21 4 ...
  .. .. .. .. .. ..@ factors : list()
  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  .. .. .. .. .. ..@ i       : int [1:5954964] 52 58 71 80 88 103 107 176 197 201 ...
  .. .. .. .. .. ..@ p       : int [1:3224] 0 1544 3756 5145 6744 9504 10322 11496 13834 17232 ...
  .. .. .. .. .. ..@ Dim     : int [1:2] 33425 3223
  .. .. .. .. .. ..@ Dimnames:List of 2
  .. .. .. .. .. .. ..$ : chr [1:33425] "DDX11L1" "WASH7P" "MIR1302-2HG" "ENSG00000238009" ...
  .. .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..@ x       : num [1:5954964] 0.86 1.63 0.86 2.06 0.86 ...
  .. .. .. .. .. ..@ factors : list()
  .. .. .. ..@ scale.data   : num [1:2000, 1:3223] -0.0345 -0.0291 -0.0569 2.8124 -0.0868 ...
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:2000] "ENSG00000235146" "LINC02593" "TNFRSF4" "HES5" ...
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. ..@ key          : chr "rna_"
  .. .. .. ..@ assay.orig   : NULL
  .. .. .. ..@ var.features : chr [1:2000] "RMST" "TENM2" "OTX2-AS1" "HTR2C" ...
  .. .. .. ..@ meta.features:'data.frame':	33425 obs. of  5 variables:
  .. .. .. .. ..$ vst.mean                 : num [1:33425] 0.00031 0.02265 0.00155 0.00434 0.00279 ...
  .. .. .. .. ..$ vst.variance             : num [1:33425] 0.00031 0.02214 0.00155 0.00433 0.00279 ...
  .. .. .. .. ..$ vst.variance.expected    : num [1:33425] 0.00031 0.02855 0.00166 0.00484 0.00304 ...
  .. .. .. .. ..$ vst.variance.standardized: num [1:33425] 1 0.776 0.936 0.894 0.916 ...
  .. .. .. .. ..$ vst.variable             : logi [1:33425] FALSE FALSE FALSE FALSE FALSE FALSE ...
  .. .. .. ..@ misc         : list()
  ..@ meta.data   :'data.frame':	3223 obs. of  4 variables:
  .. ..$ orig.ident  : Factor w/ 1 level "7316-2069": 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ nCount_RNA  : num [1:3223] 7335 4379 2479 3162 6032 ...
  .. ..$ nFeature_RNA: int [1:3223] 1544 2212 1389 1599 2760 818 1174 2338 3398 2336 ...
  .. ..$ pct_mito    : num [1:3223] 0.191 0.206 0.081 0.063 0.066 0 0.146 0.144 0.105 0.04 ...
  ..@ active.assay: chr "RNA"
  ..@ active.ident: Factor w/ 1 level "7316-2069": 1 1 1 1 1 1 1 1 1 1 ...
  .. ..- attr(*, "names")= chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  ..@ graphs      : list()
  ..@ neighbors   : list()
  ..@ reductions  :List of 7
  .. ..$ pcalognorm           :Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:30] 4.46 5.31 3.56 3.27 4.71 ...
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:30] "PClognorm_1" "PClognorm_2" "PClognorm_3" "PClognorm_4" ...
  .. .. .. ..@ feature.loadings          : num [1:2000, 1:30] 0.001694 0.002548 -0.000243 0.001501 0.000909 ...
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:2000] "RMST" "TENM2" "OTX2-AS1" "HTR2C" ...
  .. .. .. .. .. ..$ : chr [1:30] "PClognorm1" "PClognorm2" "PClognorm3" "PClognorm4" ...
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num [1:30] 9.46 5.4 5.25 5.1 4.31 ...
  .. .. .. ..@ key                       : chr "PClognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  .. ..$ umapndim20nn30lognorm:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:2] 1.63 -0.62 -3.28 -2.71 -1.08 ...
  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] -1.142 -0.329
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:2] "UMAPndim20nn30lognorm_1" "UMAPndim20nn30lognorm_2"
  .. .. .. ..@ feature.loadings          : num[0 , 0 ]
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num(0)
  .. .. .. ..@ key                       : chr "UMAPndim20nn30lognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  .. ..$ umapndim25nn30lognorm:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:2] -3.46 -2.07 1.21 1.06 -1.26 ...
  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] 1.156 0.392
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:2] "UMAPndim25nn30lognorm_1" "UMAPndim25nn30lognorm_2"
  .. .. .. ..@ feature.loadings          : num[0 , 0 ]
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num(0)
  .. .. .. ..@ key                       : chr "UMAPndim25nn30lognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  .. ..$ umapndim20nn20lognorm:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:2] -1.646 0.505 2.244 -0.168 3.539 ...
  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] 0.774 -0.297
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:2] "UMAPndim20nn20lognorm_1" "UMAPndim20nn20lognorm_2"
  .. .. .. ..@ feature.loadings          : num[0 , 0 ]
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num(0)
  .. .. .. ..@ key                       : chr "UMAPndim20nn20lognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  .. ..$ umapndim25nn20lognorm:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:2] -2.552 -0.624 -0.685 0.48 2.847 ...
  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] -3.728 -0.223
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:2] "UMAPndim25nn20lognorm_1" "UMAPndim25nn20lognorm_2"
  .. .. .. ..@ feature.loadings          : num[0 , 0 ]
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num(0)
  .. .. .. ..@ key                       : chr "UMAPndim25nn20lognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  .. ..$ umapndim20nn10lognorm:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:2] -3.998 -4.553 -2.794 0.438 -8.244 ...
  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] -0.548 -0.39
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:2] "UMAPndim20nn10lognorm_1" "UMAPndim20nn10lognorm_2"
  .. .. .. ..@ feature.loadings          : num[0 , 0 ]
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num(0)
  .. .. .. ..@ key                       : chr "UMAPndim20nn10lognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  .. ..$ umapndim25nn10lognorm:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  .. .. .. ..@ cell.embeddings           : num [1:3223, 1:2] 2.32 2.39 -1.04 -2.98 3.91 ...
  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] -0.51 0.719
  .. .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. .. ..$ : chr [1:3223] "7316-2069:AAACCCACAGTCGTTA" "7316-2069:AAACCCAGTGCCTTTC" "7316-2069:AAACCCAGTGCGCTCA" "7316-2069:AAACGAAAGCAACAGC" ...
  .. .. .. .. .. ..$ : chr [1:2] "UMAPndim25nn10lognorm_1" "UMAPndim25nn10lognorm_2"
  .. .. .. ..@ feature.loadings          : num[0 , 0 ]
  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ]
  .. .. .. ..@ assay.used                : chr "RNA"
  .. .. .. ..@ global                    : logi(0)
  .. .. .. ..@ stdev                     : num(0)
  .. .. .. ..@ key                       : chr "UMAPndim25nn10lognorm_"
  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ]
  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ]
  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ]
  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ]
  .. .. .. ..@ misc                      : list()
  ..@ images      : list()
  ..@ project.name: chr "proj"
  ..@ misc        : list()
  ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
  .. ..$ : int [1:3] 4 1 3
  ..@ commands    :List of 3
  .. ..$ NormalizeData.RNA       :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  .. .. .. ..@ name       : chr "NormalizeData.RNA"
  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2023-12-04 21:40:15"
  .. .. .. ..@ assay.used : chr "RNA"
  .. .. .. ..@ call.string: chr [1:2] "NormalizeData(seurat_obj, normalization.method = \"LogNormalize\", " "    nfeatures = 2000, assay = \"RNA\")"
  .. .. .. ..@ params     :List of 5
  .. .. .. .. ..$ assay               : chr "RNA"
  .. .. .. .. ..$ normalization.method: chr "LogNormalize"
  .. .. .. .. ..$ scale.factor        : num 10000
  .. .. .. .. ..$ margin              : num 2
  .. .. .. .. ..$ verbose             : logi TRUE
  .. ..$ FindVariableFeatures.RNA:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  .. .. .. ..@ name       : chr "FindVariableFeatures.RNA"
  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2023-12-04 21:40:18"
  .. .. .. ..@ assay.used : chr "RNA"
  .. .. .. ..@ call.string: chr [1:2] "FindVariableFeatures(seurat_obj, selection.method = \"vst\", " "    nfeatures = 2000)"
  .. .. .. ..@ params     :List of 12
  .. .. .. .. ..$ assay              : chr "RNA"
  .. .. .. .. ..$ selection.method   : chr "vst"
  .. .. .. .. ..$ loess.span         : num 0.3
  .. .. .. .. ..$ clip.max           : chr "auto"
  .. .. .. .. ..$ mean.function      :function (mat, display_progress)
  .. .. .. .. ..$ dispersion.function:function (mat, display_progress)
  .. .. .. .. ..$ num.bin            : num 20
  .. .. .. .. ..$ binning.method     : chr "equal_width"
  .. .. .. .. ..$ nfeatures          : num 2000
  .. .. .. .. ..$ mean.cutoff        : num [1:2] 0.1 8
  .. .. .. .. ..$ dispersion.cutoff  : num [1:2] 1 Inf
  .. .. .. .. ..$ verbose            : logi TRUE
  .. ..$ ScaleData.RNA           :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  .. .. .. ..@ name       : chr "ScaleData.RNA"
  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2023-12-04 21:40:21"
  .. .. .. ..@ assay.used : chr "RNA"
  .. .. .. ..@ call.string: chr "ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))"
  .. .. .. ..@ params     :List of 10
  .. .. .. .. ..$ features          : chr [1:2000] "RMST" "TENM2" "OTX2-AS1" "HTR2C" ...
  .. .. .. .. ..$ assay             : chr "RNA"
  .. .. .. .. ..$ model.use         : chr "linear"
  .. .. .. .. ..$ use.umi           : logi FALSE
  .. .. .. .. ..$ do.scale          : logi TRUE
  .. .. .. .. ..$ do.center         : logi TRUE
  .. .. .. .. ..$ scale.max         : num 10
  .. .. .. .. ..$ block.size        : num 1000
  .. .. .. .. ..$ min.cells.to.block: num 3000
  .. .. .. .. ..$ verbose           : logi TRUE
  ..@ tools       : list()
```