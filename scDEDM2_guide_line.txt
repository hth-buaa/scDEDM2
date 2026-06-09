##############################################################################################
#                                                                                            #
#    A Guide to Constructing eGRNs Using scDEDM with Paired scRNA-seq and scATAC-seq Data    #
#                                                                                            #
##############################################################################################

# Use R v4.4




########## Preparations: Installing R Packages and Setting Up the Environment ##########
### Environment Setup
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

### Environment Configuration
# Based on our development and testing, igraph version 2.0.3 and Seurat version >5.0 are required.
# While strict version requirements for other packages are currently unspecified, compatibility is most likely with versions released between 2024 and 2026.
# Nevertheless, this tutorial provides a complete list of the specific package versions used by the developers for your reference.

install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz", repos = NULL, type = "source")
packageVersion("igraph") # [1] ‘2.0.3’

install.packages(c("dplyr", "ggplot2", "Signac", "SeuratObject", "Seurat", "VGAM", "ggrepel",
                   "minpack.lm", "data.table", "tibble", "tidyr", "pROC", "DDRTree", "GA", "psych"))
packageVersion("dplyr") # [1] ‘1.1.4’
packageVersion("ggplot2") # [1] ‘4.0.0’
packageVersion("Signac") # [1] ‘1.15.0’
packageVersion("SeuratObject") # [1] ‘5.2.0’
packageVersion("Seurat") # [1] ‘5.3.0’
packageVersion("VGAM") # [1] ‘1.1.13’
packageVersion("ggrepel") # [1] ‘0.9.6’
packageVersion("minpack.lm") # [1] ‘1.2.4’
packageVersion("data.table") # [1] ‘1.17.8’
packageVersion("tibble") # [1] ‘3.3.0’
packageVersion("tidyr") # [1] ‘1.3.1’
packageVersion("pROC") # [1] ‘1.19.0.1’
packageVersion("DDRTree") # [1] ‘0.1.5’
packageVersion("GA") # [1] ‘3.2.4’
packageVersion("psych") # [1] ‘2.5.6’

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomeInfoDb", "TFBSTools", "JASPAR2024", "IRanges", "ChIPseeker",
                       "Biobase", "BiocGenerics", "monocle", "Biostrings"))
packageVersion("GenomeInfoDb") # [1] ‘1.42.3’
packageVersion("TFBSTools") # [1] ‘1.44.0’
packageVersion("JASPAR2024") # [1] ‘0.99.6’
packageVersion("IRanges") # [1] ‘2.40.1’
packageVersion("ChIPseeker") # [1] ‘1.42.1’
packageVersion("Biobase") # [1] ‘2.66.0’
packageVersion("BiocGenerics") # [1] ‘0.52.0’
packageVersion("monocle") # [1] ‘2.34.0’
packageVersion("Biostrings") # [1] ‘2.74.1’

# Load the transcript database (TxDb) package corresponding to the analysis data.
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene") # hg38-GRCh38-2013
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene") # hg19-GRCh37-2010
packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene") # [1] ‘3.20.0’
packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene") # [1] ‘3.2.2’

# Load the Ensembl genome database R annotation package based on the analysis data.
BiocManager::install("EnsDb.Hsapiens.v86") # hg38-GRCh38-2013
# BiocManager::install("EnsDb.Hsapiens.v75") # hg19-GRCh37-2010
packageVersion("EnsDb.Hsapiens.v86") # [1] ‘2.99.0’
packageVersion("EnsDb.Hsapiens.v75") # [1] ‘2.99.0’

# Load the BSgenome data package based on the analysis data.
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") # hg38-GRCh38-2013
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") # hg19-GRCh37-2010
packageVersion("BSgenome.Hsapiens.UCSC.hg38") # [1] ‘1.4.5’
packageVersion("BSgenome.Hsapiens.UCSC.hg19") # [1] ‘1.4.3’

devtools::install_github("hth-buaa/scDEDM2")
packageVersion("scDEDM") # [1] ‘2.0’

library(scDEDM)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v86)
# library(EnsDb.Hsapiens.v75)
library(dplyr)
library(DDRTree)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Hsapiens.UCSC.hg19)

### Path Setup
# The working directory for all following code execution will remain fixed at the currently set path and will not change.
Working_directory = "/mnt/i/scDEDM" # In WSL2. Strongly recommended.
# Working_directory = "I:/scDEDM" # In Windows. Not recommended.
# Working_directory = "/home/scDEDM" # In Linux. Available.
setwd(Working_directory)

### Tips
# scDEDM utilizes parallel computing, and in some cases, system memory may not be fully released after a parallel step completes.
# Restarting R can help free up memory.
# During development, each step of scDEDM was encapsulated as a function.
# To facilitate memory release, each scDEDM function can run independently—sequential execution is not enforced, but all intermediate results must be saved.
# This tutorial will demonstrate this approach and will load the necessary input data and packages before each standalone step.




########## Part 0: Loading Data ##########
##### 0.1 Loading the Seurat Object of the Example Dataset #####
setwd(Working_directory)


##### 0.1.1 Loading the Raw Seurat Object of the Example Dataset
# During development, version 0.2.2.9002 of the R package SeuratData was used.
devtools::install_github('satijalab/seurat-data')
packageVersion("SeuratData") # [1] '0.2.2.9002'

# If the download fails, try a few more times. It should succeed when the network speed is good.
options(timeout = 999999999)
SeuratData::InstallData("pbmcMultiome")
# ...
# trying URL 'http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.4.tar.gz'
# Content type 'application/octet-stream' length 280456325 bytes (267.5 MB)
# ...
scRNAseq = SeuratData::LoadData("pbmcMultiome", "pbmc.rna")
scATACseq = SeuratData::LoadData("pbmcMultiome", "pbmc.atac")

# If the R code consistently fails to download the data, you can also resort to a manual download.
# Manually download the file from the URL (http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.2.2.tar.gz) to current working directory,
# then load them.
# This tutorial uses the data obtained using this method.
install.packages("./pbmcMultiome.SeuratData_0.1.4.tar.gz", repos = NULL, type = "source")
packageVersion("pbmcMultiome.SeuratData") # [1] ‘0.1.4’
scRNAseq = pbmcMultiome.SeuratData::pbmc.rna
scATACseq = pbmcMultiome.SeuratData::pbmc.atac

base::print(scRNAseq)
# An object of class Seurat
# 36601 features across 11909 samples within 1 assay
# Active assay: RNA (36601 features, 0 variable features)
# 2 layers present: counts, data
base::print(scATACseq)
# Loading required package: Signac
# An object of class Seurat
# 108377 features across 11909 samples within 1 assay
# Active assay: ATAC (108377 features, 0 variable features)
# 2 layers present: counts, data

# Ensure that the scRNA-seq and scATAC-seq datasets contain the same set of cells.
base::print(base::table(base::sort(base::colnames(scRNAseq)) == base::sort(base::colnames(scATACseq))))
#  TRUE
# 11909

# If your scRNA-seq and scATAC-seq data come from different experiments, batches, or were generated via separate single-omics sequencing runs,
# you will need to adjust accordingly.
# For example, use appropriate methods to establish a one-to-one correspondence between cells across the two omics,
# discard cells without a matched counterpart, and then unify the cell identifiers.
# This tutorial does not provide code for unifying cell names.


##### 0.1.2 Upgrading the Seurat v4 Object to a Seurat v5 Object
# scDEDM was developed based on Seurat v5, so the format of scRNAseq needs to be converted to make it compatible with scDEDM.
if (class(scRNAseq[["RNA"]]) == "Assay") {
  scRNAseq[["RNA"]] = methods::as(object = scRNAseq[["RNA"]], Class = "Assay5")
}
base::print(class(scRNAseq[["RNA"]]))
# [1] "Assay5"
# attr(,"package")
# [1] "SeuratObject"
base::print(class(scATACseq[["ATAC"]]))
# [1] "ChromatinAssay"
# attr(,"package")
# [1] "Signac"


##### 0.1.3 Removing Data from Abnormal Chromosomes
# View chromosome names in the data.
base::print(utils::head(base::rownames(scATACseq)))
# [1] "chr1-10109-10357"   "chr1-180730-181630" "chr1-191491-191736" "chr1-267816-268196" "chr1-586028-586373" "chr1-629721-630172"
base::print(base::as.data.frame(base::table(base::sapply(base::strsplit(base::rownames(scATACseq), "-"), `[`, 1))))
#          Var1  Freq
# 1  GL000194.1     4
# 2  GL000195.1     2
# 3  GL000205.2     6
# 4  GL000219.1    10
# 5  KI270713.1     6
# 6  KI270721.1     1
# 7  KI270726.1     2
# 8  KI270727.1     1
# 9  KI270734.1     1
# 10       chr1 10538
# 11      chr10  4967
# 12      chr11  5249
# 13      chr12  5670
# 14      chr13  2410
# 15      chr14  3734
# 16      chr15  3513
# 17      chr16  4033
# 18      chr17  5766
# 19      chr18  2013
# 20      chr19  4881
# 21       chr2  8556
# 22      chr20  3153
# 23      chr21  1370
# 24      chr22  2445
# 25       chr3  6808
# 26       chr4  4268
# 27       chr5  5332
# 28       chr6  6652
# 29       chr7  5281
# 30       chr8  4416
# 31       chr9  4414
# 32       chrX  2872
# 33       chrY     3

# If abnormal chromosomes are present in the data, their associated data should be removed.
peaks_keep = Signac::seqnames(Signac::granges(scATACseq)) %in% GenomeInfoDb::standardChromosomes(Signac::granges(scATACseq))
base::print(peaks_keep)
# logical-Rle of length 108377 with 2 runs
#   Lengths: 108344     33
#   Values :   TRUE  FALSE

scATACseq = scATACseq[base::as.vector(peaks_keep), ]; rm(peaks_keep)


##### 0.1.4 Adding Cell Labels
# The cell names in the scRNA-seq and scATAC-seq data must be identical.
# Additionally, both datasets should have an identical column named seurat_annotations in their metadata
# (i.e., scRNAseq@meta.data$seurat_annotations and scATACseq@meta.data$seurat_annotations),
# representing the cell type annotations for each cell (the annotations per cell should also match).
base::print(base::as.data.frame(base::table(scRNAseq@meta.data$seurat_annotations)))
base::print(base::as.data.frame(base::table(scATACseq@meta.data$seurat_annotations)))
#              Var1 Freq
# 1       CD14 Mono 2812
# 2       CD16 Mono  514
# 3       CD4 Naive 1419
# 4         CD4 TCM 1149
# 5         CD4 TEM  298
# 6       CD8 Naive 1410
# 7       CD8 TEM_1  325
# 8       CD8 TEM_2  358
# 9            HSPC   26
# 10 Intermediate B  353
# 11           MAIT  137
# 12       Memory B  371
# 13             NK  468
# 14        Naive B  142
# 15         Plasma   18
# 16           Treg  162
# 17            cDC  198
# 18       filtered 1497
# 19            gdT  146
# 20            pDC  106

utils::head(scRNAseq@meta.data$seurat_annotations)
utils::head(scATACseq@meta.data$seurat_annotations)
# [1] "CD4 Naive" "CD4 TCM"   "CD4 Naive" "filtered"  "CD8 Naive" "CD4 Naive"

# If your data's meta.data (class: dataframe) does not contain a seurat_annotations column,
# or if the content of the seurat_annotations column does not represent cell labels,
# you will need to either manually add the cell labels to a seurat_annotations column
# or adjust the content of the existing seurat_annotations column to be the cell labels (this column should be a character vector).


##### 0.1.5 Sampling the Cells
# If your data contains multiple cell types and you are only interested in one or a few of them,
# you should filter out the cell types of no interest.

# Based on the cell labels in the data, potential cell differentiation pathways can be inferred, such as:
# 1. Myeloid Pathway
# Myeloid cells are primarily responsible for innate immunity.
# Key developmental relationships in the data include:
# (1) Monocyte Differentiation:
# CD14 Mono (classical monocytes) can differentiate into CD16 Mono (non-classical monocytes),
# which have enhanced patrolling functions.
# (2) Dendritic Cell Branch:
# cDC (conventional dendritic cells) and pDC (plasmacytoid dendritic cells) arise from a common myeloid precursor and represent distinct functional subsets.
# 2. Lymphoid Pathway
# Lymphoid cells mediate adaptive immunity, and their pathways are more clearly defined:
# (1) B-cell Lineage:
# This is a typical linear differentiation process.
# Upon antigen encounter, Naive B cells first become Intermediate B cells,
# which subsequently differentiate into either antibody-secreting Plasma cells or long-term immune Memory B cells, depending on signals received.
# (2) T-cell Lineage:
# ▪ CD4+ T-cell Path:
# Activated CD4 Naive cells can differentiate into CD4 TCM (central memory T cells),
# which may further convert to CD4 TEM (effector memory T cells).
# A separate subset of CD4 Naive cells differentiates into specialized Tregs (regulatory T cells).
# ▪ CD8+ T-cell Path:
# Activated CD8 Naive cells directly differentiate into cytotoxic effector cells, such as CD8 TEM_1/2 (effector memory CD8 T cells).

# For example, if you are only interested in the CD4 Naive cells in the dataset, you simply need to retain all CD4 Naive cells.
# Here, we choose to keep CD14 Mono and CD16 Mono for demonstration purposes.
BiocManager::install("S4Vectors")
packageVersion("S4Vectors") # [1] ‘0.44.0’
scRNAseq = S4Vectors::subset(scRNAseq, subset = seurat_annotations %in% base::c("CD14 Mono", "CD16 Mono"))
scATACseq = S4Vectors::subset(scATACseq, subset = seurat_annotations %in% base::c("CD14 Mono", "CD16 Mono"))

# Ensure that the scRNA-seq and scATAC-seq datasets contain the same set of cells.
base::print(base::table(base::sort(base::colnames(scRNAseq)) == base::sort(base::colnames(scATACseq))))
# TRUE
# 3326

# Regardless of how many cell labels of interest are selected,
# the values in the seurat_annotations column of the final data's meta.data must be unique.
# Here, we have selected two cell labels: "CD14 Mono" and "CD16 Mono", which are collectively referred to as "Mono".
scRNAseq@meta.data$original_seurat_annotations = scRNAseq@meta.data$seurat_annotations
scATACseq@meta.data$original_seurat_annotations = scATACseq@meta.data$seurat_annotations
scRNAseq@meta.data$seurat_annotations = "Mono"
scATACseq@meta.data$seurat_annotations = "Mono"


##### 0.1.6 Downsampling the Genes (Optional Step)
# To expedite runtime, downsampling can be applied to the data (selectively reducing its volume)
# to lessen the computational burden while striving to preserve the biological characteristics of the original dataset.
# However, downsampling incurs a loss of partial information, which may diminish the ability to capture subtle biological signals.

# This tutorial focuses on demonstrating the core workflow of scDEDM.
# Downsampling will be applied to cells, genes, and chromatin peaks to expedite your familiarization with the workflow.

### Cell downsampling strategy:
# first compute the non‑zero expression frequency per gene and the non‑zero chromatin accessibility frequency per peak for each cell.
# Sum these two frequencies to obtain a combined non‑zero frequency, then retain the top 3000 cells ranked by this combined frequency.
# However, we discourage excessive downsampling of cells.
# Since subsequent steps involve pseudotime trajectory analysis and grouping cells into metacells to mitigate feature sparsity,
# using too few cells would result in fewer groups and a less comprehensive gene repertoire.
# This may ultimately compromise your downstream analyses.
rna_counts = SeuratObject::GetAssayData(scRNAseq, assay = "RNA", layer = "counts")
rna_nonzero_freq = Matrix::colSums(rna_counts > 0) / base::nrow(rna_counts)
base::print(utils::head(rna_nonzero_freq))
# AAACAGCCATCCAGGT-1 AAACCAACACAATGCC-1 AAACCAACAGGAACTG-1 AAACCAACATAATCCG-1 AAACCAACATTGTGCA-1 AAACCGAAGCTGGACC-1
#         0.08272998         0.06691074         0.06650092         0.09196470         0.06576323         0.07546242

atac_counts = SeuratObject::GetAssayData(scATACseq, assay = "ATAC", layer = "counts")
atac_nonzero_freq = Matrix::colSums(atac_counts > 0) / base::nrow(atac_counts)
base::print(utils::head(atac_nonzero_freq))
# AAACAGCCATCCAGGT-1 AAACCAACACAATGCC-1 AAACCAACAGGAACTG-1 AAACCAACATAATCCG-1 AAACCAACATTGTGCA-1 AAACCGAAGCTGGACC-1
#         0.07938603         0.06072325         0.08798826         0.12506461         0.07789079         0.07502031

base::print(base::table(base::names(rna_nonzero_freq) == base::names(atac_nonzero_freq)))
# TRUE
# 3326

combined_freq = rna_nonzero_freq + atac_nonzero_freq
top_cells = base::names(combined_freq)[base::order(combined_freq, decreasing = TRUE)][1:3000]

scRNAseq = S4Vectors::subset(scRNAseq, cells = top_cells)
scATACseq = S4Vectors::subset(scATACseq, cells = top_cells)
rm(rna_counts, rna_nonzero_freq, atac_counts, atac_nonzero_freq, combined_freq, top_cells)

base::print(base::table(scRNAseq@meta.data[["original_seurat_annotations"]]))
base::print(base::table(scATACseq@meta.data[["original_seurat_annotations"]]))
# CD14 Mono CD16 Mono
#      2501       499

### Gene downsampling strategy:
# first filter out specific non-coding RNAs (those whose names begin with "AC", "AL", or "LINC") and mitochondrial genes to reduce noise.
# Then, based on the normalized data, select the top 10,000 genes with the highest expression variability.
genes = base::setdiff(base::rownames(scRNAseq), base::grep("(^AC|^AL|^LINC)[0-9]+", base::rownames(scRNAseq), value = TRUE))
genes = base::grep("MT-", genes, invert = TRUE, value = TRUE)
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)

scRNAseq = Seurat::NormalizeData(scRNAseq)
scRNAseq = Seurat::FindVariableFeatures(scRNAseq, selection.method = "vst", nfeatures = base::nrow(scRNAseq))

RNA_var_meta = scRNAseq@assays[["RNA"]]@meta.data
base::rownames(RNA_var_meta) = base::rownames(scRNAseq)
RNA_var_meta = RNA_var_meta[base::order(RNA_var_meta$vf_vst_counts_variance.standardized, decreasing = TRUE),]
genes = base::rownames(RNA_var_meta)[1:10000]
scRNAseq = S4Vectors::subset(scRNAseq, features = genes)
rm(RNA_var_meta, genes)

### Peak downsampling strategy:
# based on the variability of each peak across cells, select the top 50,000 most variable peaks.
scATACseq = Signac::RunTFIDF(scATACseq)
scATACseq = Signac::FindTopFeatures(scATACseq, min.cutoff = "q0")

ATAC_var_meta = scATACseq@assays[["ATAC"]]@meta.features
ATAC_var_meta = ATAC_var_meta[base::order(ATAC_var_meta$percentile, decreasing = TRUE), ]
peaks = base::rownames(ATAC_var_meta)[1:50000]
scATACseq = S4Vectors::subset(scATACseq, features = peaks)
rm(ATAC_var_meta, peaks)


##### 0.1.7 Saving the Example Dataset
# Use saveRDS to preserve data.
# This tutorial typically employs base::saveRDS and base::readRDS for data storage.
base::saveRDS(scRNAseq, file = "./scRNAseq.rds")
base::saveRDS(scATACseq, file = "./scATACseq.rds")

# For larger datasets, you can utilize qs::qsave to store and qs::qread to load data to reduce runtime.
install.packages("qs")
packageVersion("qs") # [1] ‘0.27.3’
qs::qsave(scRNAseq, file = "./scRNAseq.qsave")
qs::qsave(scATACseq, file = "./scATACseq.qsave")



##### 0.2 Downloading the fragment file for scATAC-seq data of the Example Dataset #####
# The fragment file for scATAC-seq data and its tbi file are essential for the Signac framework to analyze scATAC‑seq data,
# and the development of scDEDM also relies on Signac.

# If the source data provides both the _fragments.tsv.gz file and the _fragments.tsv.gz.tbi file, you can skip this step.
# However, in general, some source data (such as the data in this tutorial) only provides _fragments.tsv.gz,
# in which case we need to generate the corresponding TBI file (_fragments.tsv.gz.tbi) ourselves.
# Unfortunately, there are also many datasets that do not provide either _fragments.tsv.gz or _fragments.tsv.gz.tbi,
# in which case we have to generate both ourselves.
# This tutorial provides code to generate both _fragments.tsv.gz and _fragments.tsv.gz.tbi;
# you should choose the appropriate steps based on the specifics of your own data.

# In general, the _fragments.tsv.gz file generated by this method differs from the source data but maintains a reasonable level of similarity.

setwd(Working_directory)
scATACseq = readRDS("./scATACseq.rds")
# scATACseq = qs::qread("./scATACseq.qsave")


##### 0.2.1 Get _fragments.tsv.gz
### Download the _fragments.tsv.gz file.
# The URL of the fragment file for the scATAC‑seq data of the Example Dataset is stored in scATACseq@assays[["ATAC"]]@fragments[[1]]@path.
base::print(scATACseq@assays[["ATAC"]]@fragments[[1]]@path)
# [1] "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# Download this file to current working directory.
# This tutorial uses the the fragment file above.

# If your scATAC-seq data does not have a corresponding _fragments.tsv.gz file (for instance: the source does not provide this file, which is a common scenario.),
# use the following method to generate one.

### Generate the _fragments.tsv.gz file.
# Generate the file of fragments and its tbi file.
count_matrix = scATACseq@assays[["ATAC"]]@counts
base::print(class(count_matrix))
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"

fragments = base::rownames(count_matrix)
base::print(utils::head(fragments))
# [1] "chr19-39381822-39438328" "chr19-1230013-1281741"   "chr19-13824929-13854962" "chr14-75267507-75300328" "chr17-7832492-7859042"   "chr19-12774552-12797715"

cells = base::colnames(count_matrix)
base::print(utils::head(cells))
# [1] "AAACAGCCATCCAGGT-1" "AAACCAACACAATGCC-1" "AAACCAACAGGAACTG-1" "AAACCAACATAATCCG-1" "AAACCAACATTGTGCA-1" "AAACCGAAGCTGGACC-1"

# In this dataset, each fragment record consists of three parts separated by hyphens:
# the first is the chromosome name (formatted by concatenating "chr" with the chromosome number like 1, 2, ..., 22, X, Y),
# the second is the start position of the genomic fragment,
# and the third is the end position.
# Here, these three pieces of information are extracted separately:
# the chromosome name is stored as a character vector, while the fragment start and end positions are both stored as integer vectors.
# If the format of the fragment records in your dataset differs from this, please adjust it to match this format before proceeding with the subsequent code.
fragments_split = base::strsplit(fragments, "-")
chromosome = base::sapply(fragments_split, function(x) x[1])
start = as.integer(base::sapply(fragments_split, function(x) x[2]))
end = as.integer(base::sapply(fragments_split, function(x) x[3]))
base::print(base::as.data.frame(base::table(chromosome)))
#    chromosome Freq
# 1  GL000195.1    1
# 2  GL000205.2    2
# 3  GL000219.1    5
# 4  KI270713.1    1
# 5        chr1 4913
# 6       chr10 2260
# 7       chr11 2558
# 8       chr12 2743
# 9       chr13 1016
# 10      chr14 1604
# 11      chr15 1677
# 12      chr16 1943
# 13      chr17 2939
# 14      chr18  854
# 15      chr19 2639
# 16       chr2 3723
# 17      chr20 1506
# 18      chr21  630
# 19      chr22 1169
# 20       chr3 3059
# 21       chr4 1853
# 22       chr5 2328
# 23       chr6 2983
# 24       chr7 2369
# 25       chr8 1952
# 26       chr9 2048
# 27       chrX 1225

# Generate and save the _fragments.tsv file.
triplet = DOSE::summary(count_matrix)
result_df = data.frame(
  chromosome = chromosome[triplet$i],
  start = start[triplet$i],
  end = end[triplet$i],
  cell = cells[triplet$j],
  counts = triplet$x,
  stringsAsFactors = FALSE
)
result_df = result_df[base::order(result_df$chromosome, result_df$start, result_df$end, result_df$cell), ]
result_df = result_df[grepl("chr", result_df$chromosome), ] # Remove non‑standard chromosomes.
utils::head(result_df)
#       chromosome  start    end               cell counts
# 2977        chr1 777634 779926 AAACAGCCATCCAGGT-1      2
# 10401       chr1 777634 779926 AAACCAACACAATGCC-1      2
# 27206       chr1 777634 779926 AAACCAACATAATCCG-1      2
# 45677       chr1 777634 779926 AAACCGAAGCTGGACC-1      2
# 61061       chr1 777634 779926 AAACGCGCAGCAAGAT-1      2
# 77069       chr1 777634 779926 AAAGCACCACAATTAC-1      2

write.table(result_df, file = "./pbmc_atac_fragments.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
rm(count_matrix, fragments, cells, fragments_split, chromosome, start, end, triplet, result_df)


##### 0.2.2 Generating _fragments.tsv.gz.tbi
### Run in Linux
# Install tabix tool (tabix is a tool for indexing and fast retrieval of large data files)
sudo apt-get install tabix

# Change to working directory
cd /mnt/i/scDEDM

# Create backup of original file
cp pbmc_atac_fragments.tsv pbmc_atac_fragments_bak.tsv

# Compress TSV file to GZIP format
# bgzip is a special version of gzip that creates compressed files suitable for tabix indexing.
# The -i parameter creates an index file.
bgzip -i pbmc_atac_fragments.tsv

# Create tabix index for compressed file
# -p bed specifies the input file is in BED format.
# tabix automatically creates an index file with the same name as the compressed file but with .tbi extension.
tabix -p bed pbmc_atac_fragments.tsv.gz




########## Part 1: Data Preprocessing ##########
##### 1.1 Annotating Chromatin Fragments and Identifying TGs #####
# This step performs genomic annotation on chromatin accessibility peaks (from scATAC-seq) to identify potential target genes (TGs) with open promoter.
# It uses the ChIPseeker package for peak annotation and generates output files including annotation plots (peak_anno.png), annotation results (peak_anno.csv), and identified TGs (TGs.csv).

setwd(Working_directory)

# The genome assembly for this dataset is hg38.
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

# If the genome assembly for your dataset is hg19.
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# genome_for_anno = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

# Manually set the promoter range according to your research needs (unit: b), such as 50, 100, 200, 500, 1000, 3000, 5000, etc.
# Small range corresponds to less overall code running time, but more chromatin fragments will be ignored.
# In practice, the term "promoter" here is used in a broad sense as a collective term for both promoters and enhancers,
# i.e., cis-regulatory elements that promote regulation.
# For demonstration purposes, the promoter search range is limited to 50 base pairs upstream and downstream of the transcription start site (TSS) to reduce computational load.
# Typically, a larger range should be set, such as 3000 base pairs.
# It is important to note that subsequent steps will also involve this parameter, so its value must remain consistent throughout the entire process.
promoter_range = 50

scRNAseq = readRDS("./scRNAseq.rds")
scATACseq = readRDS("./scATACseq.rds")

library(dplyr)
results_identify_TGs = scDEDM::identify_TGs(
  genome_for_anno = genome_for_anno,
  promoter_range = promoter_range,
  scRNAseq = scRNAseq,
  scATACseq = scATACseq,
  annoDb = "org.Hs.eg.db"
)
# Time difference of 1.038718 mins

# For details on parsing the function output, refer to the function's help documentation.
?identify_TGs

base::saveRDS(results_identify_TGs, file = "./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")



##### 1.2 Identifying TFs #####
# This step extracts transcription factors Genes (TFs, TFGs) from the JASPAR2024 database that are present in the provided single-cell RNA sequencing data.
# It filters TFs based on species and collection parameters, and returns a vector of TF names available in both JASPAR2024 and the input data.

setwd(Working_directory)
scRNAseq = base::readRDS("./scRNAseq.rds")

TFs_in_JASPAR2024 = scDEDM::get_TGs_from_JASPAR2024(
  scRNAseq = scRNAseq,
  species = "Homo sapiens",
  collection = base::c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED")
)
# Time difference of 17.94591 secs

# For details on parsing the function output, refer to the function's help documentation.
?get_TGs_from_JASPAR2024

base::saveRDS(TFs_in_JASPAR2024, file = "./1.2 Data Preprocessing - Get TGs From JASPAR2024/TFs_in_JASPAR2024.rds")



##### 1.3 Generating Counts Matrices for TGs/TFs Expression and TGs Activity (Seurat Object) #####
# This function processes single-cell multi-omics data to extract expression matrices for target genes (TGs) and transcription factors (TFs),
# and computes gene activity matrices from scATAC-seq data for GRN (Gene Regulatory Network) construction.

setwd(Working_directory)

# The genome assembly for this dataset is hg38.
library(EnsDb.Hsapiens.v86)
Ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

# If the genome assembly for your dataset is hg19.
# library(EnsDb.Hsapiens.v75)
# Ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75

results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
TFs_in_JASPAR2024 = base::readRDS("./1.2 Data Preprocessing - Get TGs From JASPAR2024/TFs_in_JASPAR2024.rds")
scRNAseq = base::readRDS("./scRNAseq.rds")
scATACseq = base::readRDS("./scATACseq.rds")

# The setting of promoter_rangeshould remain consistent with the previous configuration, and this applies hereafter as well (if it is referenced again).
promoter_range = 50

# Absolute path to the fragment tbi. file (both .tsv.gz and .tbi index must exist in the same directory).
path = paste0(Working_directory, "/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz") # Use the _fragments.tsv.gz file provided with the source data.
# path = paste0(Working_directory, "/pbmc_atac_fragments.tsv.gz") # Otherwise use self-generated one instead.

Basic_Info = scDEDM::get_expression_and_activity_matrix(
  TGs = results_identify_TGs$TGs,
  TFs = TFs_in_JASPAR2024,
  scRNAseq = scRNAseq,
  scATACseq = scATACseq,
  Ensdb = Ensdb,
  promoter_range = promoter_range,
  path = path
)
# Time difference of 14.51895 mins

# For details on parsing the function output, refer to the function's help documentation.
?get_expression_and_activity_matrix

base::saveRDS(Basic_Info, file = "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.rds")




########## Part 2: Inferring Pseudotime, Branching, and Grouping Cells ##########
##### 2.1 Pseudotime Analysis and Cell Branch Partitioning #####
##### 2.1.1 Obtaining Relevant Data of Cell Types of Interest
# This step processes single-cell multi-omics data for specified cell types by performing rigorous quality control,
# filtering genes and cells based on expression thresholds,
# and ensuring consistency between RNA expression and chromatin accessibility data.

setwd(Working_directory)

Basic_Info = base::readRDS("./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.rds")
interest_cell_type = base::unique(Basic_Info[["scRNAseq_for_GRN"]]@meta.data[["seurat_annotations"]])

# Please configure the settings based on the total available RAM. During execution,
# if the RAM is exhausted, it indicates that too many cores are being used, and the number should be reduced.
install.packages("parallel")
packageVersion("parallel") # [1] ‘4.4.1’
ncores = parallel::detectCores() - 1 # In Linux. During developer runtime, parallel::detectCores() returned 79.
# ncores = 1 # In windows.

interest_cell_type_data = scDEDM::get_interest_cell_type_data(
  interest_cell_type = interest_cell_type,
  Basic_Info = Basic_Info,
  ncores = ncores
)
# Time difference of 1.82433 secs

# For details on parsing the function output, refer to the function's help documentation.
?get_interest_cell_type_data

base::saveRDS(interest_cell_type_data, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.rds")


##### 2.1.2 Performing Pseudotemporal Ordering and Branching Analysis of Cell Types of Interest Based on Transcriptomic Data
# This step performs pseudotime analysis using Monocle2 and enables interactive branch partitioning based on trajectory visualization results.
# It processes scRNA-seq data for multiple cell types to reconstruct developmental trajectories and identify branching points.

setwd(Working_directory)
interest_cell_type_data = base::readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.rds")

library(DDRTree)
interest_cell_type_Branches = scDEDM::order_pseudotime_and_divide_branches(
  interest_cell_type_data = interest_cell_type_data
)
# ...
# Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): 2
# Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
# Please enter the integer vector of cell state order for branch 1 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
# 1 2 3 9 4 5
# Please enter the integer vector of cell state order for branch 2 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
# 1 2 3 9 4 6 7 8
# ...
# Time difference of 1.706556 mins

# For details on parsing the function output, refer to the function's help documentation.
?order_pseudotime_and_divide_branches

base::saveRDS(interest_cell_type_Branches, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.rds")



##### 2.2 Cell Grouping and Group-Wise Information Specification #####
##### 2.2.1 Retrieving Pseudotemporal Expression/Activity Vectors (TFE_t, TGA_t, TGE_t) for Each Gene Based on the Pseudotime Value of Individual Cells
# This step processes single-cell multi-omics data to extract temporal expression and activity profiles for transcription factors (TFs) and target genes (TGs) along developmental trajectories.
# It performs rigorous quality filtering, maps cells to pseudotime indices, normalizes data,
# and returns organized temporal profiles for downstream GRN inference.

setwd(Working_directory)
interest_cell_type_Branches = base::readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.rds")
interest_cell_type_data = base::readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.rds")

# Manually set two gene filtering thresholds to retain genes (TFs and TGs) based on importance scores (quantile cutoffs).
# Lower values for thess parameters compress the data more aggressively, reducing subsequent runtime, but risk removing important genes.
alpha_TF = 0.69 # Numeric (0 to 1). Threshold for TF filtering based on TF Importance Score (TFIS). TFs with TFIS above this quantile will be retained. Default is 0.7.
alpha_TG = 0.62 # Numeric (0 to 1). Threshold for TG filtering based on TG Importance Score (TGIS). TGs with TGIS above this quantile will be retained. Default is 0.7.

# Manually set a cell filtering threshold to retain cells with a missing value rate below this threshold.
# Setting this parameter value too low is not recommended, as it would filter out an excessive number of cells.
# This is because scDEDM merges multiple cells into a meta‑cell and should retain informative cells rather than discard them.
# Therefore, the main purpose of this parameter is to remove overly sparse, low‑quality cells (which may not be normal biological cells).
alpha_cell = 0.98 # Numeric (0 to 1). Threshold for cell filtering based on zero-expression rate. Cells with zero-expression rate higher than this value will be removed. Default is 0.9.

# Manually set weight parameters for TF importance score (TFIS) calculation.
# Each weight ranges from -0.2 to 1, with positive values indicating positive contribution.
# Features include: non-zero expression ratio (NER), its normalized rank (NERR),
# normalized non-zero expression average (NEA), its normalized rank (NEAR),
# and the corresponding four metrics for activity data (NAR, NARR, NAA, NAAR).
w_TFNER = 0.66 # Numeric (-0.2 to 1). Weight for TF Non-Zero Expression Ratio in TFIS calculation.
w_TFNERR = 1.00 # Numeric (-0.2 to 1). Weight for TF Non-zero Expression Ratio normalized Rank (max corresponds to 1, min to 0) in TFIS calculation.
w_TFNEA = 0.64 # Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Expression Average in TFIS calculation.
w_TFNEAR = -0.10 # Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Expression Average Rank (max corresponds to 1, min to 0) in TFIS calculation.
w_TFNAR = 0.03 # Numeric (-0.2 to 1). Weight for TF Non-Zero Activity Ratio in TFIS calculation.
w_TFNARR = -0.02 # Numeric (-0.2 to 1). Weight for TF Non-zero Activity Ratio normalized Rank (max corresponds to 1, min to 0) in TFIS calculation.
w_TFNAA = 0.51 # Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Activity Average in TFIS calculation.
w_TFNAAR = 0.04 # Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Activity Average Rank (max corresponds to 1, min to 0) in TFIS calculation.

# Manually set weight parameters for TG importance score (TGIS) calculation.
# Each weight ranges from -0.2 to 1, with positive values indicating positive contribution.
# Features mirror those used for TFIS (see above).
w_TGNER = 0.24 # Numeric (-0.2 to 1). Weight for TG Non-Zero Expression Ratio in TGIS calculation.
w_TGNERR = 0.87 # Numeric (-0.2 to 1). Weight for TG Non-zero Expression Ratio normalized Rank (max corresponds to 1, min to 0) in TGIS calculation.
w_TGNEA = 0.15 # Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Expression Average in TGIS calculation.
w_TGNEAR = 0.03 # Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Expression Average Rank (max corresponds to 1, min to 0) in TGIS calculation.
w_TGNAR = 0.43 # Numeric (-0.2 to 1). Weight for TG Non-Zero Activity Ratio in TGIS calculation.
w_TGNARR = 0.82 # Numeric (-0.2 to 1). Weight for TG Non-zero Activity Ratio normalized Rank (max corresponds to 1, min to 0) in TGIS calculation.
w_TGNAA = 0.30 # Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Activity Average in TGIS calculation.
w_TGNAAR = -0.02 # Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Activity Average Rank (max corresponds to 1, min to 0) in TGIS calculation.

# In addition to manually setting these parameter values directly,
# one can also leverage auxiliary methods such as prior knowledge
# (for instance, following the approach used in the corresponding paper for the scDEDM benchmark).

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

interest_cell_type_genes_pseudotime_info = scDEDM::get_genes_pseudotime_info(
  interest_cell_type_Branches = interest_cell_type_Branches,
  interest_cell_type_data = interest_cell_type_data,
  alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
  w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
  w_TFNAR = w_TFNAR, w_TFNARR = w_TFNARR, w_TFNAA = w_TFNAA, w_TFNAAR = w_TFNAAR,
  w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
  w_TGNAR = w_TGNAR, w_TGNARR = w_TGNARR, w_TGNAA = w_TGNAA, w_TGNAAR = w_TGNAAR,
  ncores = ncores
)
# Time difference of 2.455402 secs

# For details on parsing the function output, refer to the function's help documentation.
?get_genes_pseudotime_info

base::saveRDS(interest_cell_type_genes_pseudotime_info, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")


##### 2.2.2 Grouping Cells Based on Pseudotime Values and Counts of Valid (Non-Zero) Expression and Activity Values to Calculate TFE_T, TGA_T, and TGE_T per Group
# This step performs cell grouping along pseudotime trajectories based on gene expression and activity profiles.
# It implements an adaptive grouping algorithm that ensures each group contains sufficient cells with non-zero expression/activity values for reliable downstream analysis.
# This step uses a logarithmic fitting model to determine optimal group sizes and produces aggregated expression values for each cell group.

setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")

# Manually designate three points in the first quadrant,
# where the x-axis represents the total number of cells and the y-axis represents the minimum number of cells per group.
# Use the method of undetermined coefficients to determine the parameters of the function y = c * ln(ax + b),
# fitting the relationship between the total number of cells and the minimum number of cells per group.

# Manually set a numeric vector of length 3 requiring all positive integers in strictly increasing order.
points_x_for_fitting_nc_and_nmin = base::c(200, 2500, 5000)
points_y_for_fitting_nc_and_nmin = base::c(4, 8, 12)

# If points_y_for_fitting_nc_and_nmin forms an arithmetic sequence,
# it is recommended that twice the median value of points_x_for_fitting_nc_and_nmin is less than the sum of its first and last values.
# Suppose the cell type you wish to analyze has 1,234 cells, and you want to directly set its n_min to 12.
# In this case, simply select the point (1234, 12). This means setting the value at the i-th position of points_x_for_fitting_nc_and_nmin to 1234, and the value at the i-th position of points_y_for_fitting_nc_and_nmin to 12.

ncores = parallel::detectCores() - 1 # In Linux.
# ncores = 1 # In windows

interest_cell_type_group = scDEDM::cell_grouping(
  interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
  points_x_for_fitting_nc_and_nmin = points_x_for_fitting_nc_and_nmin,
  points_y_for_fitting_nc_and_nmin = points_y_for_fitting_nc_and_nmin,
  ncores = ncores
)
# Time difference of 31.91271 secs

# For details on parsing the function output, refer to the function's help documentation.
?cell_grouping

# Check the number of cell groups in each branch, as well as the minimum non-zero counts of gene expression and gene activity within each group of cells in every branch.
# If you are unsatisfied with the grouping results
# (due to too few groups or insufficient statistical significance caused by low minimum non-zero counts of gene expression and activity),
# you need to readjust points_x_for_fitting_nc_and_nminand points_y_for_fitting_nc_and_nmin.
base::unlist(interest_cell_type_group[[1]][["Branches_n_group"]])
# [1] 7 7
base::unlist(interest_cell_type_group[[1]][["n_min"]])
# c c
# 8 8

base::saveRDS(interest_cell_type_group, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")




########## Part 3: Constructing Initial GRNs from TFBS PWMs ##########
##### 3.1 Constructing an Initial GRN for Comprehensive TG-TF Pairs Based on Base PWM of TF Binding Sites #####
# This step constructs initial gene regulatory networks (iGRNs) by predicting transcription factor binding sites (TFBS) using position weight matrices (PWMs) from the JASPAR2024 database.
# It maps transcription factors (TFs) to target genes (TGs) based on PWM matching scores in promoter regions,
# generating TF-TG association matrices for each cell type.
setwd(Working_directory)
interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")

promoter_range = 50 # Maintain consistency with earlier settings.

# The genome assembly for this dataset is hg38.
library(BSgenome.Hsapiens.UCSC.hg38)
genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# If the genome assembly for your dataset is hg19.
# library(BSgenome.Hsapiens.UCSC.hg19)
# genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# Manually set the minimum score threshold (a percentage string between 0 and 1) for PWM matching TFG.
# Although manual configuration is allowed during development, it is recommended to set it to "0%".
min_score_for_matchPWM = "0%"

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

library(dplyr)
interest_cell_type_iGRN = scDEDM::get_iGRN_by_TFBS_pwm_by_JASPAR2024(
  interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
  promoter_range = promoter_range,
  results_identify_TGs = results_identify_TGs,
  genome = genome,
  min_score_for_matchPWM = min_score_for_matchPWM,
  species = "Homo sapiens",
  collection = base::c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
  output_predicted_TFBS = FALSE,
  ncores = ncores
)
# Time difference of 1.807372 mins

# Further, users may also customize the initial regulatory strength (iGRN value) for specific TF–TG pairs based on their own research—for instance,
# by incorporating additional prior knowledge.
# (for instance, following the approach used in the corresponding paper for the scDEDM benchmark).

# For details on parsing the function output, refer to the function's help documentation.
?get_iGRN_by_TFBS_pwm_by_JASPAR2024

base::saveRDS(interest_cell_type_iGRN, file = "./3 get iGRN/interest_cell_type_iGRN.rds")



##### 3.2 Setting the Regulatory Threshold tau #####
# In the GRN, TG-TF pairs exceeding this threshold tau are considered to have a regulatory relationship;
# otherwise, they are considered to have no regulatory relationship.

# This step computes a cell-type-specific threshold (tau) for binarizing GRN.
# The threshold is calculated as the minimum non-zero interaction strength in the iGRN minus a small constant, with a lower bound of 0.005.

setwd(Working_directory)
interest_cell_type_iGRN = readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")

interest_cell_type_tau = scDEDM::get_cell_type_tau(
  interest_cell_type_iGRN = interest_cell_type_iGRN
)
# Time difference of 0.006854057 secs

base::print(interest_cell_type_tau)
#      Mono
# 0.2348684

# For details on parsing the function output, refer to the function's help documentation.
?get_cell_type_tau

base::saveRDS(interest_cell_type_tau, file = "./3 get iGRN/interest_cell_type_tau.rds")



##### 3.3 Constructing the iGRN for Each Branch #####
# This step extracts branch-specific initial gene regulatory networks (iGRNs) from the complete cell-type-specific iGRNs.
# It creates subnetworks for each branch within each cell type by filtering the global TF-TG association matrix to include only the TFs and TGs present in each specific branch.

setwd(Working_directory)
interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

interest_cell_type_branch_iGRN = scDEDM::get_branch_iGRN(
  interest_cell_type_iGRN = interest_cell_type_iGRN,
  interest_cell_type_group = interest_cell_type_group,
  ncores = ncores
)
# Time difference of 0.2698174 secs

# For details on parsing the function output, refer to the function's help documentation.
?get_branch_iGRN

base::saveRDS(interest_cell_type_branch_iGRN, file = "./3 get iGRN/interest_cell_type_branch_iGRN.rds")




########## Part 4: Developing DEDM-Based GRN Regression Predictive Models for Each Branch ##########
##### 4.1 Partitioning Training Sets #####
# This step extracts branch-specific initial gene regulatory networks (iGRNs) from the complete cell-type-specific iGRNs.
# It creates subnetworks for each branch within each cell type by filtering the global TF-TG association matrix to include only the TFs and TGs present in each specific branch.
# Substep1: Generates all possible TF–TG pairs from the iGRN and attaches the regulatory strength (theta_i).
# Substep2: For each TF and each TG, independently selects the top \code{top_n_for_max_theta_in_each_regulon} interactions (by theta_i) to form an initial training set.
# Substep3: Applies three threshold filters to the initial training set:
# (a) a quantile threshold on the thetas within the training set (\code{max_theta_retention_rate_in_tops}),
# (b) a quantile threshold on the thetas from all interactions (\code{max_theta_retention_rate_in_all}),
# and (c) an absolute number threshold (\code{retention_number_of_tops}).
# The final training set must satisfy all three.
# Substep4: If \code{refine} is TRUE,
# further refines the training set by keeping only those interactions where the TF appears in the top regulons of the TG and the TG appears in the top regulons of the TF,
# and then retains the top \code{top_ratio_in_refine} of thetas.
# Substep5: For each branch, attaches the corresponding cell group state features (TFE, TGA, TGE) for each time point to the selected interactions.
# Substep6: Returns a list for each cell type, with each element corresponding to a branch (and an "all" element) containing the training set data frame.

setwd(Working_directory)
interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")

# Manually set the number of top theta values to keep for each TF and each TG independently (a positive integer).
# Generally, it is set to 1; otherwise, the computational load becomes too large.
# If the data scale is very small, you can try setting a value greater than 1.
# top_n_for_max_theta_in_each_regulon = 1 # This is advisable.
top_n_for_max_theta_in_each_regulon = 2 # This is used for demonstration purposes in this tutorial.

# Manually set the quantile threshold applied to thetas within the initial training set (a numeric in [0, 1]).
# Here, it is set to 1, which does not reduce the default training set size.
max_theta_retention_rate_in_tops = 1

# Manually set the quantile threshold applied to thetas from all possible interactions (a numeric in [0, 1]).
# Here, it is set to 1, which does not reduce the default training set size.
max_theta_retention_rate_in_all = 1

# Manually set the number of top theta interactions to retain (a positive integer).
# Here, it is set to a large number, which does not reduce the default training set size.
retention_number_of_tops = 9999999

# Manually set this logical.
# If TRUE, further refines the training set by requiring mutual top membership (TF in top regulons of TG, and TG in top regulons of TF) and retains a top ratio.
refine = TRUE

# Manually set this logical.
# When refine = TRUE, if TRUE, ensure that each TF is included in the training set.
keep_all_TF = TRUE

# Manually set this numeric in [0, 1]. When refine = TRUE, retains the top proportion of thetas.
# Here, it is set to 1, which does not reduce the default training set size.
top_ratio_in_refine = 1

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

interest_cell_type_branch_training_set = scDEDM::select_training_set(
  interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
  interest_cell_type_group = interest_cell_type_group,
  top_n_for_max_theta_in_each_regulon = top_n_for_max_theta_in_each_regulon,
  max_theta_retention_rate_in_tops = max_theta_retention_rate_in_tops,
  max_theta_retention_rate_in_all = max_theta_retention_rate_in_all,
  retention_number_of_tops = retention_number_of_tops,
  refine = refine,
  keep_all_TF = keep_all_TF,
  top_ratio_in_refine = top_ratio_in_refine,
  ncores = ncores
)
# Time difference of 1.363743 secs

# Display the number of training sets.
lapply(interest_cell_type_branch_training_set[["Mono"]], function(x) nrow(x))
# $branch1
# [1] 93
#
# $branch2
# [1] 100
#
# $all
# [1] 105

# For details on parsing the function output, refer to the function's help documentation.
?scDEDM::select_training_set

# Of course, users can manually add training sets themselves.

# In addition, users can also manually remove entries from the training set.
# For example, if a TF in the training set corresponds to multiple TGs, and one wishes to keep only the TF–TG pair with the highest theta_i for each TF, the following code should be executed.

### Keep only highest theta_i per TF to reduce computational cost
n_branch = length(interest_cell_type_branch_training_set[[1]]) - 1
for (n in 1:n_branch) {
  df = interest_cell_type_branch_training_set[[1]][[paste0("branch", n)]]
  interest_cell_type_branch_training_set[[1]][[paste0("branch", n)]] = df[df$theta_i == ave(df$theta_i, df$TF, FUN = max), ]
}

# Combine unique TF-TG pairs across all branches
combined = base::unique(base::do.call(base::rbind, base::lapply(interest_cell_type_branch_training_set[[1]][-n_branch - 1], function(x) x[, 1:3])))
base::rownames(combined) = paste0(combined$TF, "_to_", combined$TG)
interest_cell_type_branch_training_set[[1]][["all"]] = combined
rm(n_branch, n, combined)

# Display the number of training sets.
lapply(interest_cell_type_branch_training_set[["Mono"]], function(x) nrow(x))
# $branch1
# [1] 39
#
# $branch2
# [1] 40
#
# $all
# [1] 43

# This tutorial uses this additionally filtered training set for subsequent demonstrations.

base::saveRDS(interest_cell_type_branch_training_set, file = "./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")



##### 4.2 Build Predictive Model for Each Branch #####
##### 4.2.1 Setting Initial Values, Lower Bounds, and Upper Bounds for Parameters
# This step initializes the parameter structures for the gene regulatory network prediction model scDEDM.

setwd(Working_directory)
interest_cell_type_branch_training_set = base::readRDS("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")

interest_cell_type_branch_init_params = scDEDM::set_init_params(
  interest_cell_type_branch_training_set = interest_cell_type_branch_training_set
)
# Time difference of 0.7103155 secs

# For details on parsing the function output, refer to the function's help documentation.
?scDEDM::set_init_params

base::saveRDS(interest_cell_type_branch_init_params, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds")


##### 4.2.2 Training models
# This step trains prediction models for gene regulatory networks using a hybrid optimization approach combining genetic algorithms and gradient ascent.
# It processes each branch independently, training individual models for each regulator-target pair (regulon).
# The training incorporates multiple early stopping criteria and parallel processing capabilities for efficient computation.

setwd(Working_directory)
interest_cell_type_branch_training_set = base::readRDS("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")
interest_cell_type_branch_init_params = base::readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds")
interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")

### Manually set maximum training epochs.
# This parameter is generally set to a relatively large value,
# as training termination is controlled by an early‑stopping strategy rather than by reaching the maximum number of training rounds.
# It is advisable to first set this parameter to a moderate value (such as 50) to run the model and observe the training epochs required for each training set.
# In some cases, it is acceptable for a small number of training sets to reach the maximum training epochs.
# On one hand, the results from these training sets are not necessarily poor; on the other hand, this approach helps save overall runtime.
# Directly setting a very high number of training epochs often leads to significantly longer total training time due to a few difficult-to-train sets,
# but the cost-effectiveness of such extra time investment tends to be low.
max_epochs = 50

### Manually set early stopping epochs.
# A higher value for this parameter leads to longer runtime, which may result in better training performance but also increases the risk of overfitting.
early_stop = 5

### Manually set parameters related to the genetic algorithm.
# Probability of crossover between chromosome pairs in genetic algorithm.
# It is not recommended to modify this parameter.
# This is because the training process is primarily driven by gradient descent,
# while the genetic algorithm is used to help the training escape local optima (serving as a corrective mechanism).
# Setting this parameter in this way corresponds to a stronger corrective force,
# which, although it sacrifices some of the genetic algorithm’s contribution to training progress, is negligible,
# given that the convergence speed of gradient descent far exceeds that of the genetic algorithm.
pcrossover_ga = 0.5

# Probability of mutation in parent chromosomes in genetic algorithm.
# It is not recommended to modify this parameter.
# The reasoning is consistent with the explanation above.
pmutation_ga = 0.5

# Population size for each genetic algorithm training session.
# A higher value for this parameter increases the runtime,
# while also improving its effectiveness in escaping local optima.
popSize_ga = 50

# Number of generations for the genetic algorithm per iteration round.
# A higher value for this parameter increases the runtime,
# while also improving its effectiveness in escaping local optima.
maxiter_ga = 150

### Manually set parameters related to the gradient descent method.
# Number of gradient descent steps per iteration round.
# A higher value for this parameter results in longer runtime,
# while the marginal effect on training progress diminishes accordingly.
iterations_grad = 30

### Manually set parameters related to result selection.
# Threshold for filtering results based on difference between predicted and observed regulatory strengths.
# A lower value for this parameter increases the runtime (as it imposes stricter requirements and makes convergence more difficult).
theta_difference_threshold = 0.1

# Threshold for filtering results based on fitting loss.
# A lower value for this parameter increases the runtime (as it imposes stricter requirements and makes convergence more difficult).
fit_loss_threshold = 0.1

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

interest_cell_type_branch_model_train = scDEDM::model_train(
  interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
  interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
  interest_cell_type_group = interest_cell_type_group,
  interest_cell_type_tau = interest_cell_type_tau,
  interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
  max_epochs = max_epochs, early_stop = early_stop,
  eps_theta = 1e-3, eps_loss = 1e-5,
  popSize_ga = popSize_ga, maxiter_ga = maxiter_ga, pcrossover_ga = pcrossover_ga, pmutation_ga = pmutation_ga,  parallel_ga = FALSE, seed_ga = 123,
  iterations_grad = iterations_grad, ncores_grad = 1,
  theta_difference_threshold = theta_difference_threshold, fit_loss_threshold = fit_loss_threshold,
  ncores = ncores
)
# ...
# Starting training the prediction model for branch 1 in cell type Mono.
# There are 39 regulons for training.
# ...
# Starting training the prediction model for branch 2 in cell type Mono.
# There are 40 regulons for training.
# ...
# Time difference of 9.557524 mins

# For details on parsing the function output, refer to the function's help documentation.
?model_train

base::saveRDS(interest_cell_type_branch_model_train, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")




########## Part 5: Predicting GRN ##########
##### 5.1 Getting Data for Inferring GRN #####
# This step prepares comprehensive datasets for gene regulatory network (GRN) inference by combining iGRN with cell group state features.
# It generates all possible TF-TG pairs, adds their corresponding regulatory strength values,
# and integrates cell group state measurements (TFE, TGA, TGE) for each time point.
# Data is processed in parallel for efficiency and can be partitioned to handle large datasets.

setwd(Working_directory)
interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

# Manually set the number of partitions to slice the data for parallel processing.
# It is recommended to partition the data into as many subsets as the number of cores utilized during execution.
n_part = ncores

interest_cell_type_branch_data_for_GRN_infer = scDEDM::get_data_for_inferring_GRN(
  interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
  interest_cell_type_group = interest_cell_type_group,
  n_part = ncores,
  ncores = ncores
)
# ...
# Process branch 1 in cell type Mono.
# There are 1394 potential regulatory relationship.
# ...
# Process branch 2 in cell type Mono.
# There are 1394 potential regulatory relationship.
# ...
# Time difference of 2.42798 secs

# For details on parsing the function output, refer to the function's help documentation.
?get_data_for_inferring_GRN

base::saveRDS(interest_cell_type_branch_data_for_GRN_infer, file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")



##### 5.2 Inferring GRN for Each Branch #####
# This step infers gene regulatory networks (GRN) by applying trained models to predict regulatory strengths for all potential TF-TG interactions.
# It processes data in parallel, calculates predicted regulatory strengths (theta_p) and fitting losses, and organizes results by cell type and branch.
# It supports selective processing of specific cell types and branches for targeted analysis.

setwd(Working_directory)
interest_cell_type_branch_data_for_GRN_infer = base::readRDS("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")
interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_branch_model_train = base::readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

### Way 1: one step (If the running memory is sufficient)
# Manually set the maximum number of regulatory interactions to process in each parallel chunk.
# It is recommended to distribute the data load as evenly as possible across each core.
each = ceiling((max(base::sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][1:(length(interest_cell_type_branch_data_for_GRN_infer[[1]]) - 1)],
  base::nrow
))) / ncores)

interest_cell_type_pGRN = scDEDM::infer_GRN(
  interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
  interest_cell_type_tau = interest_cell_type_tau,
  interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
  interest_cell_type_group = interest_cell_type_group,
  interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
  select_cell_type = NULL,
  select_branch = NULL,
  each = each,
  ncores = ncores
)
# ...
# Process branch 1 in cell type Mono.
# 1394 potential regulatory relationships will be predicted.
# 39 trained regulatory relationships will be predicted.
# ...
# Process branch 2 in cell type Mono.
# 1394 potential regulatory relationships will be predicted.
# 40 trained regulatory relationships will be predicted.
# ...
# Time difference of 1.141221 mins

# If the dataset is large or contains many branches, this method can easily run out of memory.
# In such cases, it is advisable to use Way 2 described below.
# For example, when running the benchmark dataset and the apply dataset with scDEDM, this method can cause memory explosion.
# Therefore, the Way 2 approach is adopted, where each branch is processed separately.

### Way 2: more steps but stable
# If the process is terminated due to memory issues,
# you can manually continue running it branch by branch (saving the results to interest_cell_type_pGRN).
interest_cell_type_pGRN = list()
cell_type = base::names(interest_cell_type_branch_data_for_GRN_infer)[1]
interest_cell_type_pGRN[[cell_type]] = list()

# If executing the entire for-loop code still risks a memory explosion, or if you are concerned about this possibility,
# manually set the parameter n (n=1, n=2, ...) and run the code separately for each branch.
# The reason for this issue is that with large or multi-branch datasets,
# system memory may not be fully released promptly after each parallel computing task,
# leading to a gradual accumulation of memory usage.
for (n in 1:(length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
  message("Infer cell type ", cell_type, " branch ", n, ".")
  each = ceiling(base::nrow(interest_cell_type_branch_data_for_GRN_infer[[1]][[paste0("branch", n)]]) / ncores)
  interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]] = scDEDM::infer_GRN(
    interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tau = interest_cell_type_tau,
    interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
    interest_cell_type_group = interest_cell_type_group,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    select_cell_type = cell_type,
    select_branch = n,
    each = each,
    ncores = ncores
  )[[cell_type]][[paste0("branch", n)]]
}
# branch 1 - 1394 potential regulatory relationships - 39 trained regulatory relationships - 33.43238 secs
# branch 2 - 1394 potential regulatory relationships - 40 trained regulatory relationships - 32.47773 secs

# For details on parsing the function output, refer to the function's help documentation.
?infer_GRN

base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")



##### 5.3 Selecting Best Prediction for each TF-TG #####
# This step filters and selects the best predictions from inferred gene regulatory networks (GRN) by applying multiple quality thresholds.
# It processes each branch independently, selects the optimal regulatory strength (theta_p) and fitting loss for each TF-TG pair,
# and integrates results across branches using weighted averaging based on cell counts.

setwd(Working_directory)
interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
interest_cell_type_branch_model_train = base::readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
interest_cell_type_tau = base::readRDS("./3 get iGRN/interest_cell_type_tau.rds")
interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
interest_cell_type_branch_data_for_GRN_infer = base::readRDS("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")

# Manually set parameters related to quality control.
# You can control the network scale by adjusting these parameters.
fit_loss_quantile_threshold1 = 0.2 # Quantile threshold for filtering predictions based on loss in initial selection.
theta_p_quantile_threshold1 = 0.8 # Quantile threshold for filtering predictions based on regulatory strength in initial selection.
fit_loss_upper1 = 0.2 # Upper bound for acceptable fitting loss in initial selection.
theta_p_lower1 = 0.2 # Lower bound for acceptable regulatory strength in initial selection.

ncores = parallel::detectCores() - 1 # linux or wsl2
# cores = 1 # win

# Manually set the maximum number of regulatory interactions to process in each parallel chunk.
# It is recommended to distribute the data load as evenly as possible across each core.
each = ceiling((max(base::sapply(
  interest_cell_type_branch_data_for_GRN_infer[[1]][1:(length(interest_cell_type_branch_data_for_GRN_infer[[1]]) - 1)],
  base::nrow
))) / ncores)

interest_cell_type_pGRN = scDEDM::select_best_prediction(
  interest_cell_type_pGRN = interest_cell_type_pGRN,
  interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
  interest_cell_type_tau = interest_cell_type_tau,
  interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
  fit_loss_quantile_threshold = fit_loss_quantile_threshold1,
  theta_p_quantile_threshold = theta_p_quantile_threshold1,
  fit_loss_upper = fit_loss_upper1,
  theta_p_lower = theta_p_lower1,
  each = each,
  ncores = ncores
)
# Time difference of 7.680629 secs

# For details on parsing the function output, refer to the function's help documentation.
?select_best_prediction

base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")

# Now, we have obtained the per-branch gene regulatory networks (TF → TG) without quality control.
# The first column, theta_p, represents the predicted regulatory strength, ranging from 0 to 1.
# The second column, loss, indicates the corresponding prediction error or fitting loss.
# The third column, fit, specifies the training samples used for the fitting process.

base::print(utils::head(interest_cell_type_pGRN[[1]][["branch1"]][["best_result"]]))
#                   theta_p_branch1 loss_branch1    fit_branch1
# SREBF2_to_FOXJ3         0.4816801   0.13288868 BCL6_to_PLAGL1
# SREBF2_to_SDCCAG8       0.8553212   0.11397457 BCL6_to_PLAGL1
# SREBF2_to_SVIL          0.5751898   0.07473404 ZBTB11_to_SVIL
# SREBF2_to_GLUD1         0.6986417   0.09009428 BCL6_to_PLAGL1
# SREBF2_to_OSBPL5        0.3270159   0.14519183   PBX3_to_RHOG
# SREBF2_to_RHOG          0.1920971   0.16368286 ZBTB11_to_SVIL

base::print(utils::head(interest_cell_type_pGRN[[1]][["branch2"]][["best_result"]]))
#                   theta_p_branch2 loss_branch2          fit_branch2
# SREBF2_to_FOXJ3         0.7864455   0.23922368 SREBF2_to_CSGALNACT1
# SREBF2_to_SDCCAG8       0.6715592   0.33383910     FOXK1_to_LRRFIP2
# SREBF2_to_PIP4K2A       0.3547520   0.60741796       ZBTB11_to_SVIL
# SREBF2_to_SVIL          0.9249903   0.09386839         XBP1_to_SVIL
# SREBF2_to_GLUD1         0.3071636   0.31118698       NFIL3_to_GLUD1
# SREBF2_to_OSBPL5        0.4277840   0.07834429      NR2C2_to_OSBPL5

# Furthermore, we have obtained the quality‑controlled, branch‑aggregated gene regulatory network (TF → TG).
# Columns n‑1 through 3×n contain the GRN data from branch n.
# The second‑to‑last column, theta_p, represents the integrated predicted regulatory strength obtained by weighted summation.
# The last column, fit_loss, corresponds to the integrated prediction error or fitting loss derived from the same weighted summation for this integrated theta_p.

base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary"]]))
#                theta_p_branch1 loss_branch1      fit_branch1 theta_p_branch2 loss_branch2       fit_branch2    theta_p  fit_loss
# ATF3_to_ATP10A       0.2383470    0.2948810     XBP1_to_SVIL      0.16656969    0.1202470  TRPS1_to_ST3GAL5 0.20402820 0.2113835
# ATF3_to_BBX          0.6055109    0.3951221 SMAD3_to_LRRFIP2      0.32133603    0.3639518     RFX3_to_PSMB9 0.46963880 0.3802187
# ATF3_to_CABIN1       0.4786749    0.4018531   CEBPD_to_PIPOX      0.36453017    0.4287394    CEBPD_to_PIPOX 0.42409906 0.4147082
# ATF3_to_CCDC57       0.4895349    0.3214108   CEBPD_to_PIPOX      0.26828992    0.6882797    CEBPD_to_PIPOX 0.38375140 0.4968212
# ATF3_to_CD44                NA           NA             <NA>      0.02763036    0.3835237      XBP1_to_SVIL 0.02763036 0.3835237
# ATF3_to_COX4I1              NA           NA             <NA>      0.32572531    0.3157606 NR4A2_to_SERPINA1 0.32572531 0.3157606



##### 5.4 Filtering GRN #####
# This step applies quality thresholds to filter predicted gene regulatory networks (GRN) at two stringency levels: loose and strict.
# It uses quantile-based thresholds and absolute bounds to select high-confidence regulatory interactions based on fitting loss and predicted regulatory strength (theta_p).
# Results are saved in separate data frames for each filtering level.

setwd(Working_directory)
interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")

# Manually set parameters related to quality control.
# You can control the network scale by adjusting these parameters.
fit_loss_quantile_threshold2 = 0.3 # Quantile threshold for filtering predictions based on loss in further selection.
theta_p_quantile_threshold2 = 0.7 # Quantile threshold for filtering predictions based on regulatory strength in further selection.
fit_loss_upper2 = 0.4 # Upper bound for acceptable fitting loss in further selection.
theta_p_lower2 = 0.4 # Lower bound for acceptable regulatory strength in further selection.

interest_cell_type_pGRN = scDEDM::filter_GRN(
  interest_cell_type_pGRN = interest_cell_type_pGRN,
  fit_loss_quantile_threshold = fit_loss_quantile_threshold2,
  theta_p_quantile_threshold = theta_p_quantile_threshold2,
  fit_loss_upper = fit_loss_upper2,
  theta_p_lower = theta_p_lower2
)
# Time difference of 0.005716801 secs

# For details on parsing the function output, refer to the function's help documentation.
?filter_GRN

base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")

# We now have two gene regulatory networks (TF → TG) with different levels of quality‑control stringency.
# The "strict" version retains only regulatory relationships that are present in every branch,
# whereas the "loose" version includes all relationships found in any branch.

base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary_loose"]]))
#                 theta_p_branch1 loss_branch1     fit_branch1 theta_p_branch2 loss_branch2     fit_branch2   theta_p   fit_loss
# ATF3_to_PSMB9         0.6205020   0.05448074   ATF3_to_PSMB9       0.6198214   0.08747811   ATF3_to_PSMB9 0.6201766 0.07025772
# ATF7_to_ATP10A        0.7610145   0.17437031 CTCF_to_MAP3K14              NA           NA            <NA> 0.7610145 0.17437031
# ATF7_to_MAP3K14       0.5887166   0.08787508 CTCF_to_MAP3K14       0.9156424   0.20717070 CTCF_to_MAP3K14 0.7450291 0.14491370
# ATF7_to_PSMB9         0.6272701   0.02304706   ATF7_to_PSMB9       0.6272671   0.07924958   ATF7_to_PSMB9 0.6272687 0.04991908
# ATF7_to_SVIL          0.3388653   0.11030983  ZBTB11_to_SVIL       0.8465623   0.21871372    XBP1_to_SVIL 0.5816096 0.16214081
# BACH2_to_PSMB9               NA           NA            <NA>       0.7323872   0.03308331  BACH2_to_PSMB9 0.7323872 0.03308331

base::print(utils::head(interest_cell_type_pGRN[[1]][["pred_summary_strict"]]))
#                      theta_p_branch1 loss_branch1     fit_branch1 theta_p_branch2 loss_branch2          fit_branch2   theta_p   fit_loss
# ATF3_to_PSMB9              0.6205020   0.05448074   ATF3_to_PSMB9       0.6198214   0.08747811        ATF3_to_PSMB9 0.6201766 0.07025772
# ATF7_to_MAP3K14            0.5887166   0.08787508 CTCF_to_MAP3K14       0.9156424   0.20717070      CTCF_to_MAP3K14 0.7450291 0.14491370
# ATF7_to_PSMB9              0.6272701   0.02304706   ATF7_to_PSMB9       0.6272671   0.07924958        ATF7_to_PSMB9 0.6272687 0.04991908
# ATF7_to_SVIL               0.3388653   0.11030983  ZBTB11_to_SVIL       0.8465623   0.21871372         XBP1_to_SVIL 0.5816096 0.16214081
# BCL11A_to_CSGALNACT1       0.4935648   0.06027581 BCL11A_to_PPARG       0.4693906   0.18863870 SREBF2_to_CSGALNACT1 0.4820064 0.12164975
# BCL11A_to_MAP3K14          0.4278057   0.14496531  ZNF148_to_TJP2       0.7302362   0.12246965      CTCF_to_MAP3K14 0.5724063 0.13420950

# In the scDEDM paper, for all datasets run through scDEDM, subsequent analyses used the "loose" GRN because it retains more comprehensive results.
# The "strict" GRN is also valuable, as it retains the most core regulatory relationships—those that must be present in every branch.



##### 5.5 Inferring peak-TG and TF-peak #####
##### 5.5.1 Inferring peak-TG relationships
# This step infers regulatory relationships between chromatin accessibility peaks and target genes (TGs)
# by integrating scRNA-seq and scATAC-seq data with pre-computed gene regulatory network (GRN) predictions.
# It processes results at two stringency levels (loose and strict) to identify target genes present in the filtered GRN,
# then retrieves promoter-annotated peaks associated with these genes and calculates multiple importance scores for each peak-TG pair.

setwd(Working_directory)
scRNAseq = readRDS("./scRNAseq.rds")
scATACseq = readRDS("./scATACseq.rds")
interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")

# Manually set the weighting of individual components when calculating the importance score for each peak, each TG, and each peak-TG association.
# Based on benchmarking experience, the following configuration is recommended; however, it is fully customizable according to your preference.
weight_Accessibility_percentile = 2.5 # Weight for the accessibility percentile in importance score calculations
weight_Expression_percentile = 2.5 # Weight for the expression percentile in importance score calculations.
weight_Accessibility_non_zero_ratio_percentile = 1 # Weight for the non-zero ratio percentile of accessibility in importance score calculations.
weight_Expression_non_zero_ratio_percentile = 1 # Weight for the non-zero ratio percentile of expression in importance score calculations.

library(dplyr)
interest_cell_type_peak_TG_pred = scDEDM::infer_peak_to_TG(
  scRNAseq = scRNAseq,
  scATACseq = scATACseq,
  interest_cell_type_pGRN = interest_cell_type_pGRN,
  results_identify_TGs = results_identify_TGs,
  weight_Accessibility_percentile = weight_Accessibility_percentile,
  weight_Expression_percentile = weight_Expression_percentile,
  weight_Accessibility_non_zero_ratio_percentile = weight_Accessibility_non_zero_ratio_percentile,
  weight_Expression_non_zero_ratio_percentile = weight_Expression_non_zero_ratio_percentile
)
# Time difference of 12.63694 secs

# For details on parsing the function output, refer to the function's help documentation.
?infer_peak_to_TG

base::saveRDS(interest_cell_type_peak_TG_pred, file = "./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")

# Now, we have obtained the initial predictions of peak–TG associations for both the "strict" and "loose" networks.

base::print(utils::head(interest_cell_type_peak_TG_pred[[1]][["loose"]][["peak_TG"]]))
# # A tibble: 6 × 20
#   peak                 SYMBOL Accessibility Expression Accessibility_percen…¹ Expression_percentile Accessibility_non_ze…² Expression_non_zero_…³ Accessibility_non_ze…⁴
#   <chr>                <chr>          <dbl>      <dbl>                  <dbl>                 <dbl>                  <dbl>                  <dbl>                  <dbl>
# 1 chr15-25783511-2578… ATP10A         2.04        2.03                 0.774                 0.613                  0.121                   0.542                 0.290
# 2 chr22-24174874-2417… CABIN1         2.03        1.81                 0.742                 0.0645                 0.116                   0.430                 0.258
# 3 chr16-85798643-8580… COX4I1         0.741       2.13                 0.226                 0.742                  0.264                   1.28                  0.806
# 4 chr8-19642027-19642… CSGAL…         2.44        2.17                 0.871                 0.839                  0.0935                  0.591                 0.0968
# 5 chr1-42314420-42315… FOXJ3          2.47        2.04                 0.935                 0.645                  0.0855                  1.02                  0.0323
# 6 chr10-87092154-8709… GLUD1          0.477       2.05                 0.0645                0.677                  0.310                   0.913                 1
# # ℹ abbreviated names: ¹​Accessibility_percentile, ²​Accessibility_non_zero_ratio, ³​Expression_non_zero_ratio, ⁴​Accessibility_non_zero_ratio_percentile
# # ℹ 11 more variables: Expression_non_zero_ratio_percentile <dbl>, peak_importance_score <dbl>, peak_importance_score_percentile <dbl>, TG_importance_score <dbl>,
# #   TG_importance_score_percentile <dbl>, peak_TG_importance_score <dbl>, peak_TG_importance_score_percentile <dbl>, TG_count <int>, peak_count <int>,
# #   TG_count_percentile <dbl>, peak_count_percentile <dbl>

base::print(utils::head(interest_cell_type_peak_TG_pred[[1]][["strict"]][["peak_TG"]]))
# # A tibble: 6 × 20
#   peak                 SYMBOL Accessibility Expression Accessibility_percen…¹ Expression_percentile Accessibility_non_ze…² Expression_non_zero_…³ Accessibility_non_ze…⁴
#   <chr>                <chr>          <dbl>      <dbl>                  <dbl>                 <dbl>                  <dbl>                  <dbl>                  <dbl>
# 1 chr15-25783511-2578… ATP10A         2.04        2.03                 0.731                 0.577                  0.121                   0.542                 0.346
# 2 chr22-24174874-2417… CABIN1         2.03        1.81                 0.692                 0.0769                 0.116                   0.430                 0.308
# 3 chr8-19642027-19642… CSGAL…         2.44        2.17                 0.846                 0.808                  0.0935                  0.591                 0.115
# 4 chr1-42314420-42315… FOXJ3          2.47        2.04                 0.923                 0.615                  0.0855                  1.02                  0.0385
# 5 chr10-87092154-8709… GLUD1          0.477       2.05                 0.0769                0.654                  0.310                   0.913                 1
# 6 chr5-55850868-55851… IL31RA         2.18        1.91                 0.769                 0.231                  0.101                   0.258                 0.231
# # ℹ abbreviated names: ¹​Accessibility_percentile, ²​Accessibility_non_zero_ratio, ³​Expression_non_zero_ratio, ⁴​Accessibility_non_zero_ratio_percentile
# # ℹ 11 more variables: Expression_non_zero_ratio_percentile <dbl>, peak_importance_score <dbl>, peak_importance_score_percentile <dbl>, TG_importance_score <dbl>,
# #   TG_importance_score_percentile <dbl>, peak_TG_importance_score <dbl>, peak_TG_importance_score_percentile <dbl>, TG_count <int>, peak_count <int>,
# #   TG_count_percentile <dbl>, peak_count_percentile <dbl>

# In the presentation above, the first column represents the peak,
# the second column corresponds to the TG associated with that peak,
# and all subsequent columns describe features related to this peak‑TG pair.

base::print(base::colnames(interest_cell_type_peak_TG_pred[[1]][["loose"]][["peak_TG"]]))
base::print(base::colnames(interest_cell_type_peak_TG_pred[[1]][["strict"]][["peak_TG"]]))
#  [1] "peak"                                    "SYMBOL"                                  "Accessibility"
#  [4] "Expression"                              "Accessibility_percentile"                "Expression_percentile"
#  [7] "Accessibility_non_zero_ratio"            "Expression_non_zero_ratio"               "Accessibility_non_zero_ratio_percentile"
# [10] "Expression_non_zero_ratio_percentile"    "peak_importance_score"                   "peak_importance_score_percentile"
# [13] "TG_importance_score"                     "TG_importance_score_percentile"          "peak_TG_importance_score"
# [16] "peak_TG_importance_score_percentile"     "TG_count"                                "peak_count"
# [19] "TG_count_percentile"                     "peak_count_percentile"

# You can select peak–TG predictions further based on the features according to your analytical focus.
# For example, in the benchmark section of the scDEDM paper,
# given that all feature values are numeric and all features are positive indicators (a larger value indicates greater importance of the peak–TG pair),
# the developers assigned a weight to each feature:
# very small positive weights, zero weights, or even minimal negative weights were given to features of lesser focus,
# while positive weights were assigned to features of greater interest, with the weight value increasing in proportion to the level of emphasis.
# performed a weighted sum as the comprehensive importance score across all features to obtain a composite importance score for each peak–TG association,
# and then set a truncation threshold to filter the results.
# This tutorial now provides example code for this method, allowing you to set parameters (the feature weights and truncation threshold) as needed.
# In the paper, these parameters were set based on the gold standard data used in the benchmark analysis.

interest_cell_type_peak_TG_pred_filtered = interest_cell_type_peak_TG_pred

# Set weights for each feature by yourself.
# As an example, the peak_TG_importance_score is given the highest priority,
# followed by the peak_importance_score,
# and then the TG_importance_score,
# with the main focus placed on these three features.
weight_Accessibility_1 = 0.20
weight_Expression_1 = 0.15
weight_Accessibility_percentile_1 = 0.16
weight_Expression_percentile_1 = 0.12
weight_Accessibility_non_zero_ratio_1 = 0.16
weight_Expression_non_zero_ratio_1 = 0.12
weight_Accessibility_non_zero_ratio_percentile_1 = 0.16
weight_Expression_non_zero_ratio_percentile_1 = 0.12
weight_peak_importance_score_1 = 0.5
weight_peak_importance_score_percentile_1 = 0.60
weight_TG_importance_score_1 = 0.30
weight_TG_importance_score_percentile_1 = 0.40
weight_peak_TG_importance_score_1 = 0.70
weight_peak_TG_importance_score_percentile_1 = 0.80
weight_TG_count_1 = 0.02
weight_peak_count_1 = 0.02
weight_TG_count_percentile_1 = 0.08
weight_peak_count_percentile_1 = 0.08

# Combine all weights into a vector.
weights_1 = base::c(
  weight_Accessibility_1, weight_Expression_1,
  weight_Accessibility_percentile_1, weight_Expression_percentile_1,
  weight_Accessibility_non_zero_ratio_1, weight_Expression_non_zero_ratio_1,
  weight_Accessibility_non_zero_ratio_percentile_1, weight_Expression_non_zero_ratio_percentile_1,
  weight_peak_importance_score_1, weight_peak_importance_score_percentile_1,
  weight_TG_importance_score_1, weight_TG_importance_score_percentile_1,
  weight_peak_TG_importance_score_1, weight_peak_TG_importance_score_percentile_1,
  weight_TG_count_1, weight_peak_count_1,
  weight_TG_count_percentile_1, weight_peak_count_percentile_1
)

# Process both "loose" and "strict" filtering levels.
for (level in c("loose", "strict")) {
  # Get peak-TG information for current level.
  peak_TG_info = interest_cell_type_peak_TG_pred_filtered[[1]][[level]][["peak_TG"]]

  # Calculate composite importance score using weighted sum.
  peak_TG_info$peak_TG_IS = as.numeric(base::as.matrix(peak_TG_info[, base::c(
    "Accessibility", "Expression",
    "Accessibility_percentile", "Expression_percentile",
    "Accessibility_non_zero_ratio", "Expression_non_zero_ratio",
    "Accessibility_non_zero_ratio_percentile", "Expression_non_zero_ratio_percentile",
    "peak_importance_score", "peak_importance_score_percentile",
    "TG_importance_score", "TG_importance_score_percentile",
    "peak_TG_importance_score", "peak_TG_importance_score_percentile",
    "TG_count", "peak_count",
    "TG_count_percentile", "peak_count_percentile"
  )]) %*% weights_1)

  # Set truncation threshold to filter the results by yourself.
  # As an example, set this threshold to the 60th percentile of all peak_TG_IS.
  thrshold_peak_TG_IS = as.numeric(stats::quantile(peak_TG_info$peak_TG_IS, probs = 0.6, na.rm = TRUE))

  # Filter peak-TG pairs above threshold
  peak_TG_info_filtered = peak_TG_info[peak_TG_info$peak_TG_IS > thrshold_peak_TG_IS, ]
  TGs_filtered = base::sort(base::unique(peak_TG_info_filtered$SYMBOL))

  # Update results
  interest_cell_type_peak_TG_pred_filtered[[1]][[level]][["TGs"]] = TGs_filtered
  interest_cell_type_peak_TG_pred_filtered[[1]][[level]][["peak_TG"]] = peak_TG_info_filtered
}
rm(
  weight_Accessibility_1, weight_Expression_1,
  weight_Accessibility_percentile_1, weight_Expression_percentile_1,
  weight_Accessibility_non_zero_ratio_1, weight_Expression_non_zero_ratio_1,
  weight_Accessibility_non_zero_ratio_percentile_1, weight_Expression_non_zero_ratio_percentile_1,
  weight_peak_importance_score_1, weight_peak_importance_score_percentile_1,
  weight_TG_importance_score_1, weight_TG_importance_score_percentile_1,
  weight_peak_TG_importance_score_1, weight_peak_TG_importance_score_percentile_1,
  weight_TG_count_1, weight_peak_count_1,
  weight_TG_count_percentile_1, weight_peak_count_percentile_1,
  weights_1, peak_TG_info, thrshold_peak_TG_IS, peak_TG_info_filtered, TGs_filtered, level
)

base::saveRDS(interest_cell_type_peak_TG_pred_filtered, file = "./5 Infer GRN/interest_cell_type_peak_TG_pred_filtered.rds")


##### 5.5.2 Inferring TF-peak relationships
# This step infers regulatory relationships between transcription factors (TFs) and chromatin accessibility peaks by integrating filtered GRN predictions with peak-TG associations.
# It processes results at two stringency levels (loose and strict) to identify TFs present in predicted GRNs,
# maps them to their target genes,
# and identifies peaks associated with these target genes through promoter annotations.

setwd(Working_directory)
interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
interest_cell_type_peak_TG_pred = base::readRDS("./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")
interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")

ncores = parallel::detectCores() - 1 # in Linux
# ncores = 1 # in Windows

interest_cell_type_TF_peak_pred = scDEDM::infer_TF_to_peak(
  interest_cell_type_pGRN = interest_cell_type_pGRN,
  interest_cell_type_peak_TG_pred = interest_cell_type_peak_TG_pred,
  interest_cell_type_iGRN = interest_cell_type_iGRN,
  ncores = ncores
)
# Time difference of 1.386666 secs

# For details on parsing the function output, refer to the function's help documentation.
?infer_TF_to_peak

base::saveRDS(interest_cell_type_TF_peak_pred, file = "./5 Infer GRN/interest_cell_type_TF_peak_pred.rds")

# A transcription factor (TF), together with the peaks and target genes (TGs) it regulates, forms a regulon.
# Now, we have obtained the regulons for all transcription factors in the quality control level of "strict" and "loose".

base::print(utils::head(interest_cell_type_TF_peak_pred[[1]][["loose"]][["TF_peaks"]][["CTCF"]][["peaks"]]))
#                       peak     SYMBOL Accessibility Expression Accessibility_percentile Expression_percentile Accessibility_non_zero_ratio Expression_non_zero_ratio Accessibility_non_zero_ratio_percentile Expression_non_zero_ratio_percentile peak_importance_score peak_importance_score_percentile
# 1  chr15-25783511-25784351     ATP10A     2.0392602   2.026126                0.7741935             0.6129032                   0.12099611                 0.5423263                              0.29032258                            0.3225806              2.225806                        0.8064516
# 4   chr8-19642027-19642670 CSGALNACT1     2.4389707   2.171052                0.8709677             0.8387097                   0.09349388                 0.5912498                              0.09677419                            0.3870968              2.274194                        0.8387097
# 7   chr5-55850868-55851503     IL31RA     2.1773403   1.907809                0.8064516             0.2903226                   0.10088343                 0.2581902                              0.19354839                            0.1290323              2.209677                        0.7741935
# 8   chr3-37174623-37177523    LRRFIP2     0.7830701   2.181423                0.2580645             0.8709677                   0.25789110                 1.2637710                              0.77419355                            0.8064516              1.419355                        0.2580645
# 9  chr17-45284065-45284937    MAP3K14     2.5436536   1.950516                1.0000000             0.3870968                   0.08563634                 0.5923068                              0.06451613                            0.4193548              2.564516                        0.9677419
# 12   chr11-3100202-3102226     OSBPL5     1.4430224   1.947632                0.5161290             0.3548387                   0.18182082                 0.6199962                              0.58064516                            0.4516129              1.870968                        0.5483871
#    TG_importance_score TG_importance_score_percentile peak_TG_importance_score peak_TG_importance_score_percentile TG_count peak_count TG_count_percentile peak_count_percentile   theta_i theta_i_percentile   theta_p theta_p_percentile   fit_loss fit_loss_percentile
# 1            1.8548387                      0.5806452                 4.080645                           0.6129032        1          1                   1             0.4918033 0.7992338          0.1111111 0.9773799          1.0000000 0.15937538           0.7058824
# 4            2.4838710                      0.7419355                 4.758065                           0.8709677        1          1                   1             0.4918033 0.8117728          0.2222222 0.4893827          0.3333333 0.14454403           0.5882353
# 7            0.8548387                      0.2580645                 3.064516                           0.3225806        1          1                   1             0.4918033 0.9074112          0.7777778 0.4834694          0.2222222 0.16813183           0.8235294
# 8            2.9838710                      0.8709677                 4.403226                           0.8064516        1          1                   1             0.4918033 0.8366419          0.5555556 0.5267105          0.4444444 0.06882151           0.2352941
# 9            1.3870968                      0.3870968                 3.951613                           0.5483871        1          1                   1             0.4918033 0.9523810          1.0000000 0.9525596          0.8888889 0.02766003           0.1176471
# 12           1.3387097                      0.3225806                 3.209677                           0.3870968        1          1                   1             0.4918033 0.9389902          0.8888889 0.5935439          0.7777778 0.09734987           0.3529412

base::print(utils::head(interest_cell_type_TF_peak_pred[[1]][["strict"]][["TF_peaks"]][["CTCF"]][["peaks"]]))
#                      peak     SYMBOL Accessibility Expression Accessibility_percentile Expression_percentile Accessibility_non_zero_ratio Expression_non_zero_ratio Accessibility_non_zero_ratio_percentile Expression_non_zero_ratio_percentile peak_importance_score peak_importance_score_percentile
# 1 chr15-25783511-25784351     ATP10A     2.0392602   2.026126                0.7307692             0.5769231                   0.12099611                 0.5423263                              0.34615385                           0.26923077              2.173077                        0.7692308
# 3  chr8-19642027-19642670 CSGALNACT1     2.4389707   2.171052                0.8461538             0.8076923                   0.09349388                 0.5912498                              0.11538462                           0.34615385              2.230769                        0.8076923
# 6  chr5-55850868-55851503     IL31RA     2.1773403   1.907809                0.7692308             0.2307692                   0.10088343                 0.2581902                              0.23076923                           0.07692308              2.153846                        0.7307692
# 7  chr3-37174623-37177523    LRRFIP2     0.7830701   2.181423                0.2692308             0.8461538                   0.25789110                 1.2637710                              0.76923077                           0.80769231              1.442308                        0.2692308
# 8 chr17-45284065-45284937    MAP3K14     2.5436536   1.950516                1.0000000             0.3461538                   0.08563634                 0.5923068                              0.07692308                           0.38461538              2.576923                        0.9615385
# 9   chr11-3100202-3102226     OSBPL5     1.4430224   1.947632                0.4615385             0.3076923                   0.18182082                 0.6199962                              0.61538462                           0.42307692              1.769231                        0.5000000
#   TG_importance_score TG_importance_score_percentile peak_TG_importance_score peak_TG_importance_score_percentile TG_count peak_count TG_count_percentile peak_count_percentile   theta_i theta_i_percentile   theta_p theta_p_percentile   fit_loss fit_loss_percentile
# 1           1.7115385                      0.5384615                 3.884615                           0.5769231        1          1                   1             0.4901961 0.7992338          0.1111111 0.9773799          1.0000000 0.15937538           0.7058824
# 3           2.3653846                      0.7307692                 4.596154                           0.8076923        1          1                   1             0.4901961 0.8117728          0.2222222 0.4893827          0.3333333 0.14454403           0.5882353
# 6           0.6538462                      0.2307692                 2.807692                           0.2307692        1          1                   1             0.4901961 0.9074112          0.7777778 0.4834694          0.2222222 0.16813183           0.8235294
# 7           2.9230769                      0.8461538                 4.365385                           0.7692308        1          1                   1             0.4901961 0.8366419          0.5555556 0.5267105          0.4444444 0.06882151           0.2352941
# 8           1.2500000                      0.3076923                 3.826923                           0.5000000        1          1                   1             0.4901961 0.9523810          1.0000000 0.9525596          0.8888889 0.02766003           0.1176471
# 9           1.1923077                      0.2692308                 2.961538                           0.2692308        1          1                   1             0.4901961 0.9389902          0.8888889 0.5935439          0.7777778 0.09734987           0.3529412

# In the presentation above, the first column represents the peak associated with the TF (the TF above is CTCF),
# the second column corresponds to the TG associated with that peak,
# and all subsequent columns describe features related to this peak‑TG pair.

base::print(base::colnames(interest_cell_type_TF_peak_pred[[1]][["strict"]][["TF_peaks"]][["CTCF"]][["peaks"]]))
base::print(base::colnames(interest_cell_type_TF_peak_pred[[1]][["strict"]][["TF_peaks"]][["CTCF"]][["peaks"]]))
#  [1] "peak"                                    "SYMBOL"                                  "Accessibility"
#  [4] "Expression"                              "Accessibility_percentile"                "Expression_percentile"
#  [7] "Accessibility_non_zero_ratio"            "Expression_non_zero_ratio"               "Accessibility_non_zero_ratio_percentile"
# [10] "Expression_non_zero_ratio_percentile"    "peak_importance_score"                   "peak_importance_score_percentile"
# [13] "TG_importance_score"                     "TG_importance_score_percentile"          "peak_TG_importance_score"
# [16] "peak_TG_importance_score_percentile"     "TG_count"                                "peak_count"
# [19] "TG_count_percentile"                     "peak_count_percentile"                   "theta_i"
# [22] "theta_i_percentile"                      "theta_p"                                 "theta_p_percentile"
# [25] "fit_loss"                                "fit_loss_percentile"

# You can select TF-peak predictions further based on the features according to your analytical focus.
# For example, in the benchmark section of the scDEDM paper,
# given that all feature values are numeric，
# and that all features except those related to fit_loss are positive indicators (where a larger value indicates greater importance of the TF‑peak pair),
# while the two fit_loss‑related features are negative indicators (where a smaller value indicates greater importance of the TF‑peak pair),
# the developers assigned a weight to each positive feature:
# very small positive weights, zero weights, or even minimal negative weights were given to features of lesser focus,
# while positive weights were assigned to features of greater interest, with the weight value increasing in proportion to the level of emphasis.
# the developers assigned a negative weight to each negative feature, with the weight value decreasing in proportion to the level of emphasis.
# Then, the developers performed a weighted sum as the comprehensive importance score across all features to obtain a composite importance score for each TF-peak association,
# and then set a truncation threshold to filter the results.
# All TF‑peak pairs within every regulon share the above parameters (weights and truncation threshold).
# If all TF‑peak predictions for a particular regulon are filtered out due to an excessive truncation threshold,
# the prediction with the highest comprehensive importance score will be forcibly retained.
# This tutorial now provides example code for this method, allowing you to set the parameters as needed.
# In the paper, these parameters were set based on the gold standard data used in the benchmark analysis.

interest_cell_type_TF_peak_pred_filtered = interest_cell_type_TF_peak_pred

# Define weights for each feature by yourself
# As an example, the highest priority is given to theta_p,
# followed by fit_loss, then theta_i, and finally peak_TG_importance_score,
# with the main focus placed on these four features.
weight_Accessibility_2 = 0.06
weight_Expression_2 = 0.06
weight_Accessibility_percentile_2 = 0.04
weight_Expression_percentile_2 = 0.04
weight_Accessibility_non_zero_ratio_2 = 0.1
weight_Expression_non_zero_ratio_2 = 0.1
weight_Accessibility_non_zero_ratio_percentile_2 = 0.08
weight_Expression_non_zero_ratio_percentile_2 = 0.08
weight_peak_importance_score_2 = 0.3
weight_peak_importance_score_percentile_2 = 0.2
weight_TG_importance_score_2 = 0.2
weight_TG_importance_score_percentile_2 = 0.1
weight_peak_TG_importance_score_2 = 0.4
weight_peak_TG_importance_score_percentile_2 = 0.3
weight_TG_count_2 = 0.06
weight_peak_count_2 = 0.06
weight_TG_count_percentile_2 = 0.04
weight_peak_count_percentile_2 = 0.04
weight_theta_i_2 = 0.6
weight_theta_i_percentile_2 = 0.5
weight_theta_p_2 = 1
weight_theta_p_percentile_2 = 0.9
weight_fit_loss_2 = -0.8
weight_fit_loss_percentile_2 = -0.7

# Combine all weights into a vector
weights_2 = base::c(
  weight_Accessibility_2, weight_Expression_2,
  weight_Accessibility_percentile_2, weight_Expression_percentile_2,
  weight_Accessibility_non_zero_ratio_2, weight_Expression_non_zero_ratio_2,
  weight_Accessibility_non_zero_ratio_percentile_2, weight_Expression_non_zero_ratio_percentile_2,
  weight_peak_importance_score_2, weight_peak_importance_score_percentile_2,
  weight_TG_importance_score_2, weight_TG_importance_score_percentile_2,
  weight_peak_TG_importance_score_2, weight_peak_TG_importance_score_percentile_2,
  weight_TG_count_2, weight_peak_count_2,
  weight_TG_count_percentile_2, weight_peak_count_percentile_2,
  weight_theta_i_2, weight_theta_i_percentile_2,
  weight_theta_p_2, weight_theta_p_percentile_2,
  weight_fit_loss_2, weight_fit_loss_percentile_2
)

# Process both "loose" and "strict" filtering levels
for (level in c("loose", "strict")) {
  TFs = interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TFs"]]

  # Calculate composite importance score for each TF-peak pair
  for (tf in TFs) {
    TF_peak_info = interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]]

    # Weighted sum of all features
    TF_peak_info$TF_peak_IS = as.numeric(base::as.matrix(TF_peak_info[, base::c(
      "Accessibility", "Expression",
      "Accessibility_percentile", "Expression_percentile",
      "Accessibility_non_zero_ratio", "Expression_non_zero_ratio",
      "Accessibility_non_zero_ratio_percentile", "Expression_non_zero_ratio_percentile",
      "peak_importance_score", "peak_importance_score_percentile",
      "TG_importance_score", "TG_importance_score_percentile",
      "peak_TG_importance_score", "peak_TG_importance_score_percentile",
      "TG_count", "peak_count",
      "TG_count_percentile", "peak_count_percentile",
      "theta_i", "theta_i_percentile",
      "theta_p", "theta_p_percentile",
      "fit_loss","fit_loss_percentile"
    )]) %*% weights_2)

    interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]] = TF_peak_info
  }

  # As an example, set this threshold to the 60th percentile of all TF_peak_IS.
  all_TF_peak_IS = base::unlist(base::lapply(TFs, function(tf) {
    interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]][["TF_peak_IS"]]
  }))
  thrshold_TF_peak_IS = as.numeric(stats::quantile(all_TF_peak_IS, probs = 0.6, na.rm = TRUE))

  # Filter peaks above threshold, keep at least one peak per TF
  for (tf in TFs) {
    TF_peak_info = interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]]

    # If threshold exceeds all scores, keep the peak with highest score
    max_TF_peak_IS = max(TF_peak_info$TF_peak_IS)
    if (thrshold_TF_peak_IS > max_TF_peak_IS) {
      TF_peak_info = TF_peak_info[TF_peak_info$TF_peak_IS == max_TF_peak_IS, ]
    } else {
      TF_peak_info = TF_peak_info[TF_peak_info$TF_peak_IS >= thrshold_TF_peak_IS, ]
    }

    interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]] = TF_peak_info
  }

  # Update TF-TG relationship data frame
  tf_tg_df = base::unique(base::as.data.frame(data.table::rbindlist(base::lapply(TFs, function(tf) {data.frame(
    TF = tf,
    TG = interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]][["SYMBOL"]],
    stringsAsFactors = FALSE
  )}))))
  rownames(tf_tg_df) = paste0(tf_tg_df$TF, "_to_", tf_tg_df$TG)
  interest_cell_type_TF_peak_pred_filtered[[1]][[level]][["TF_TGs"]] = tf_tg_df
}
rm(
  weight_Accessibility_2, weight_Expression_2,
  weight_Accessibility_percentile_2, weight_Expression_percentile_2,
  weight_Accessibility_non_zero_ratio_2, weight_Expression_non_zero_ratio_2,
  weight_Accessibility_non_zero_ratio_percentile_2, weight_Expression_non_zero_ratio_percentile_2,
  weight_peak_importance_score_2, weight_peak_importance_score_percentile_2,
  weight_TG_importance_score_2, weight_TG_importance_score_percentile_2,
  weight_peak_TG_importance_score_2, weight_peak_TG_importance_score_percentile_2,
  weight_TG_count_2, weight_peak_count_2,
  weight_TG_count_percentile_2, weight_peak_count_percentile_2,
  weight_theta_i_2, weight_theta_i_percentile_2,
  weight_theta_p_2, weight_theta_p_percentile_2,
  weight_fit_loss_2, weight_fit_loss_percentile_2,
  weights_2, level, TFs, tf, TF_peak_info, all_TF_peak_IS, thrshold_TF_peak_IS, max_TF_peak_IS, tf_tg_df
)

base::saveRDS(interest_cell_type_TF_peak_pred_filtered, file = "./5 Infer GRN/interest_cell_type_TF_peak_pred_filtered.rds")



##### 5.6 Inferring TF-peak-TG eGRN #####
# This step integrates TF-peak and peak-TG predictions to construct a comprehensive enhancer-gene regulatory network (eGRN) linking transcription factors, regulatory peaks, and target genes.
# This step is not packaged as a function because different quality control (filtering) approaches may require different code implementations.
# Following the logic from the previous subsection,
# eGRN is extracted here from interest_cell_type_TF_peak_pred and interest_cell_type_peak_TG_pred (or from interest_cell_type_TF_peak_pred_filtered and interest_cell_type_peak_TG_pred_filtered).

setwd(Working_directory)

# Load filtered TF-peak and peak-TG predictions
# You can load either the prediction results filtered by weights and truncation thresholds, or the unfiltered prediction results.
TF_peak_pred = base::readRDS("./5 Infer GRN/interest_cell_type_TF_peak_pred_filtered.rds")
peak_TG_pred = base::readRDS("./5 Infer GRN/interest_cell_type_peak_TG_pred_filtered.rds")

# TF_peak_pred = base::readRDS("./5 Infer GRN/interest_cell_type_TF_peak_pred.rds")
# peak_TG_pred = base::readRDS("./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")

# Initialize structure for integrated eGRN
interest_cell_type_peGRN = list()
interest_cell_type_peGRN[[base::names(TF_peak_pred)[1]]] = list()

# Process both filtering levels: "loose" and "strict"
for (level in c("loose", "strict")) {
  # Combine all TF-peak predictions into a single data frame
  TF_peak_df = base::do.call(
    base::rbind,
    base::lapply(TF_peak_pred[[1]][[level]][["TFs"]], function(tf) {
      base::cbind(TF = tf, TF_peak_pred[[1]][[level]][["TF_peaks"]][[tf]][["peaks"]])
    })
  )

  # Get peak-TG predictions for current level
  peak_TG_df = peak_TG_pred[[1]][[level]][["peak_TG"]]

  # Align column names between TF-peak and peak-TG data frames
  TF_peak_df[, base::setdiff(base::colnames(peak_TG_df), base::colnames(TF_peak_df))] = NA
  peak_TG_df[, base::setdiff(base::colnames(TF_peak_df), base::colnames(peak_TG_df))] = NA
  peak_TG_df = peak_TG_df[, base::colnames(TF_peak_df)]

  # Create unique identifiers for peak-TG pairs
  TF_peak_df$peak_TG = paste0(TF_peak_df$peak, "_to_", TF_peak_df$SYMBOL)
  peak_TG_df$peak_TG = paste0(peak_TG_df$peak, "_to_", peak_TG_df$SYMBOL)


  # Initialize eGRN with TF-peak pairs
  eGRN = TF_peak_df
  columns = base::colnames(eGRN)

  # Merge peak-TG information into TF-peak pairs
  for (i in 1:base::nrow(eGRN)) {
    peak_tg = eGRN$peak_TG[i]
    if (peak_tg %in% peak_TG_df$peak_TG) {
      NAcolumns = columns[base::is.na(eGRN[i, ])]
      eGRN[i, NAcolumns] = peak_TG_df[peak_TG_df$peak_TG == peak_tg, NAcolumns, drop = TRUE]
    }
  }

  # Add peak-TG pairs that don't have corresponding TF-peak predictions
  for (i in 1:base::nrow(peak_TG_df)) {
    peak_tg = peak_TG_df$peak_TG[i]
    if (!(peak_tg %in% eGRN$peak_TG)) {
      eGRN = base::rbind(eGRN, peak_TG_df[i, , drop = FALSE])
    }
  }

  # Standardize column naming
  if("SYMBOL" %in% base::colnames(eGRN)) {
    base::colnames(eGRN)[base::colnames(eGRN) == "SYMBOL"] = "TG"
  }

  # Annotate regulation types
  eGRN$Significant_TF_peak = !base::is.na(eGRN$TF)
  eGRN$Significant_peak_TG = eGRN$peak_TG %in% peak_TG_df$peak_TG
  eGRN$Regulation_type = base::ifelse(
    eGRN$Significant_TF_peak & eGRN$Significant_peak_TG, "TF-peak-TG", base::ifelse(
      eGRN$Significant_TF_peak & !eGRN$Significant_peak_TG, "TF-peak", base::ifelse(
        !eGRN$Significant_TF_peak & eGRN$Significant_peak_TG, "peak-TG", NA
      )
    )
  )
  eGRN$Regulation_type = factor(eGRN$Regulation_type, levels = c("TF-peak-TG", "TF-peak", "peak-TG"), exclude = NULL)

  # Remove temporary identifier column
  eGRN$peak_TG = NULL

  # Create simplified version with key columns
  eGRN_simplify = eGRN[, base::c("TF", "peak", "TG", "Significant_TF_peak", "Significant_peak_TG", "Regulation_type")]
  eGRN_simplify$TG[eGRN_simplify$Regulation_type == "TF-peak"] = NA
  eGRN_simplify$TF[eGRN_simplify$Regulation_type == "peak-TG"] = NA

  # Store results
  interest_cell_type_peGRN[[1]][[level]][["eGRN"]] = eGRN
  interest_cell_type_peGRN[[1]][[level]][["eGRN_simplify"]] = eGRN_simplify
}
rm(TF_peak_pred, peak_TG_pred, level, TF_peak_df, peak_TG_df, columns, eGRN, i, peak_tg, NAcolumns, eGRN_simplify)

base::saveRDS(interest_cell_type_peGRN, file = "./5 Infer GRN/interest_cell_type_peGRN.rds")

# Now, we have obtained the final predicted loose eGRN and strict eGRN.
# In interest_cell_type_peGRN[[1]][[level]][["eGRN"]], each row corresponds to a regulatory relationship.
# The first column, TF, represents the transcription factor in this relationship;
# the second column, peak, denotes the associated peak;
# and the third column, TG, indicates the target gene.
# Subsequent columns describe features of this regulatory relationship.
# The third-to-last and second-to-last columns specify whether the TF‑peak component and the peak‑TG component, respectively, are strongly supported (passed previous quality control).
# The final column indicates the type of strongly supported component(s) in the relationship, which can be TF‑peak‑TG, TF‑peak, or peak‑TG.

utils::head(tibble::as_tibble(interest_cell_type_peGRN[[1]][["loose"]][["eGRN"]]))
# # A tibble: 6 × 32
#   TF     peak           TG    Accessibility Expression Accessibility_percen…¹ Expression_percentile Accessibility_non_ze…² Expression_non_zero_…³ Accessibility_non_ze…⁴
#   <chr>  <chr>          <chr>         <dbl>      <dbl>                  <dbl>                 <dbl>                  <dbl>                  <dbl>                  <dbl>
# 1 ATF3   chr6-32849879… PSMB9         0.507       2.17                 0.0968                 0.806                  0.301                  1.33                   0.935
# 2 ATF7   chr15-2578351… ATP1…         2.04        2.03                 0.774                  0.613                  0.121                  0.542                  0.290
# 3 ATF7   chr6-32849879… PSMB9         0.507       2.17                 0.0968                 0.806                  0.301                  1.33                   0.935
# 4 ATF7   chr10-2963498… SVIL          1.54        2.25                 0.548                  0.935                  0.167                  1.27                   0.452
# 5 BACH2  chr6-32849879… PSMB9         0.507       2.17                 0.0968                 0.806                  0.301                  1.33                   0.935
# 6 BCL11A chr14-9438732… SERP…         1.60        2.48                 0.613                  0.968                  0.156                  1.79                   0.419
# # ℹ abbreviated names: ¹​Accessibility_percentile, ²​Accessibility_non_zero_ratio, ³​Expression_non_zero_ratio, ⁴​Accessibility_non_zero_ratio_percentile
# # ℹ 22 more variables: Expression_non_zero_ratio_percentile <dbl>, peak_importance_score <dbl>, peak_importance_score_percentile <dbl>, TG_importance_score <dbl>,
# #   TG_importance_score_percentile <dbl>, peak_TG_importance_score <dbl>, peak_TG_importance_score_percentile <dbl>, TG_count <int>, peak_count <int>,
# #   TG_count_percentile <dbl>, peak_count_percentile <dbl>, theta_i <dbl>, theta_i_percentile <dbl>, theta_p <dbl>, theta_p_percentile <dbl>, fit_loss <dbl>,
# #   fit_loss_percentile <dbl>, TF_peak_IS <dbl>, peak_TG_IS <dbl>, Significant_TF_peak <lgl>, Significant_peak_TG <lgl>, Regulation_type <fct>

utils::head(tibble::as_tibble(interest_cell_type_peGRN[[1]][["strict"]][["eGRN"]]))
# # A tibble: 6 × 32
#   TF     peak           TG    Accessibility Expression Accessibility_percen…¹ Expression_percentile Accessibility_non_ze…² Expression_non_zero_…³ Accessibility_non_ze…⁴
#   <chr>  <chr>          <chr>         <dbl>      <dbl>                  <dbl>                 <dbl>                  <dbl>                  <dbl>                  <dbl>
# 1 ATF3   chr6-32849879… PSMB9         0.507       2.17                  0.115                 0.769                 0.301                   1.33                  0.923
# 2 ATF7   chr17-4528406… MAP3…         2.54        1.95                  1                     0.346                 0.0856                  0.592                 0.0769
# 3 ATF7   chr6-32849879… PSMB9         0.507       2.17                  0.115                 0.769                 0.301                   1.33                  0.923
# 4 BCL11A chr14-9438732… SERP…         1.60        2.48                  0.577                 0.962                 0.156                   1.79                  0.462
# 5 BCL11A chr10-2963498… SVIL          1.54        2.25                  0.5                   0.923                 0.167                   1.27                  0.5
# 6 BCL6   chr1-42314420… FOXJ3         2.47        2.04                  0.923                 0.615                 0.0855                  1.02                  0.0385
# # ℹ abbreviated names: ¹​Accessibility_percentile, ²​Accessibility_non_zero_ratio, ³​Expression_non_zero_ratio, ⁴​Accessibility_non_zero_ratio_percentile
# # ℹ 22 more variables: Expression_non_zero_ratio_percentile <dbl>, peak_importance_score <dbl>, peak_importance_score_percentile <dbl>, TG_importance_score <dbl>,
# #   TG_importance_score_percentile <dbl>, peak_TG_importance_score <dbl>, peak_TG_importance_score_percentile <dbl>, TG_count <int>, peak_count <int>,
# #   TG_count_percentile <dbl>, peak_count_percentile <dbl>, theta_i <dbl>, theta_i_percentile <dbl>, theta_p <dbl>, theta_p_percentile <dbl>, fit_loss <dbl>,
# #   fit_loss_percentile <dbl>, TF_peak_IS <dbl>, peak_TG_IS <dbl>, Significant_TF_peak <lgl>, Significant_peak_TG <lgl>, Regulation_type <fct>

base::table(interest_cell_type_peGRN[[1]][["loose"]][["eGRN"]][["Regulation_type"]])
# TF-peak-TG    TF-peak    peak-TG
#         63         19          0

base::table(interest_cell_type_peGRN[[1]][["strict"]][["eGRN"]][["Regulation_type"]])
# TF-peak-TG    TF-peak    peak-TG
#         46         19          0

# A simplified version of the eGRN is stored in interest_cell_type_peGRN[[1]][[level]][["eGRN_simply"]].

utils::head(interest_cell_type_peGRN[[1]][["loose"]][["eGRN_simplify"]])
#       TF                    peak       TG Significant_TF_peak Significant_peak_TG Regulation_type
# 1   ATF3  chr6-32849879-32856058     <NA>                TRUE               FALSE         TF-peak
# 2   ATF7 chr15-25783511-25784351   ATP10A                TRUE                TRUE      TF-peak-TG
# 3   ATF7  chr6-32849879-32856058     <NA>                TRUE               FALSE         TF-peak
# 4   ATF7 chr10-29634985-29637034     SVIL                TRUE                TRUE      TF-peak-TG
# 5  BACH2  chr6-32849879-32856058     <NA>                TRUE               FALSE         TF-peak
# 6 BCL11A chr14-94387323-94388656 SERPINA1                TRUE                TRUE      TF-peak-TG

utils::head(interest_cell_type_peGRN[[1]][["strict"]][["eGRN_simplify"]])
#       TF                    peak       TG Significant_TF_peak Significant_peak_TG Regulation_type
# 1   ATF3  chr6-32849879-32856058     <NA>                TRUE               FALSE         TF-peak
# 2   ATF7 chr17-45284065-45284937  MAP3K14                TRUE                TRUE      TF-peak-TG
# 3   ATF7  chr6-32849879-32856058     <NA>                TRUE               FALSE         TF-peak
# 4 BCL11A chr14-94387323-94388656 SERPINA1                TRUE                TRUE      TF-peak-TG
# 5 BCL11A chr10-29634985-29637034     SVIL                TRUE                TRUE      TF-peak-TG
# 6   BCL6  chr1-42314420-42315427    FOXJ3                TRUE                TRUE      TF-peak-TG
