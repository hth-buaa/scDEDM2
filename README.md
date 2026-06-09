# scDEDM: A Discrete Evolutionary Dynamical Model for eGRN Inference from Single-Cell Matched Multi-omics Data

## Project Overview
scDEDM is an R package designed to construct enhancer-gene regulatory networks (eGRNs) by integrating paired single-cell RNA sequencing (scRNA-seq) and single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq) data. To infer high-confidence regulatory relationships between transcription factors (TFs), chromatin accessibility peaks, and target genes (TGs), it leverages pseudotemporal ordering, transcription factor binding site (TFBS) prediction, a discrete evolutionary dynamical model (based on Transcription Factor Expression (TFE), Target Gene Activity (TGA), and Target Gene Expression (TGE))，and a hybrid optimization approach (combining genetic algorithms and gradient ascent). 

Key features of scDEDM include:
- Integration of multi-omics data to link chromatin accessibility with gene expression
- Pseudotime-based cell grouping to handle sparse single-cell data
- Initial network construction using TFBS position weight matrices (PWMs) from JASPAR2024
- Robust model training with early stopping to avoid overfitting
- Generation of both "loose" (comprehensive) and "strict" (core) eGRNs
- Inference of TF-peak, peak-TG, and TF-peak-TG regulatory relationships

## Installation

### Prerequisites
- **R v4.4**: scDEDM is developed and tested under R version 4.4. Ensure you have this version installed.
- **Operating System**: Compatible with Linux, Windows, and WSL2 (Windows Subsystem for Linux). Linux is recommended for parallel computing efficiency. In contrast, the use of Windows is not recommended as it does not support parallel computation.

### Installing Dependencies
First, install the required R packages. Note that igraph version 2.0.3 must be installed from the archive source.  

```r
# Install igraph v2.0.3 (required version)
install.packages("https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz", repos = NULL, type = "source")

# Install other CRAN packages
install.packages(c("dplyr", "ggplot2", "Seurat", "Signac", "VGAM", "ggrepel", "minpack.lm",
                   "data.table", "tibble", "tidyr", "pROC", "DDRTree", "GA", "psych",
                   "qs", "parallel", "devtools"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomeInfoDb", "TFBSTools", "JASPAR2024", "IRanges", "ChIPseeker",
                       "Biobase", "BiocGenerics", "monocle", "Biostrings", "DOSE"))

# Install genome annotation packages (choose one set based on your reference genome)
# For hg38/GRCh38:
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "EnsDb.Hsapiens.v86",
                       "BSgenome.Hsapiens.UCSC.hg38"))
# For hg19/GRCh37 (optional):
# BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene",
#                        "EnsDb.Hsapiens.v75",
#                        "BSgenome.Hsapiens.UCSC.hg19"))
```

### Installing scDEDM
Install scDEDM from GitHub using `devtools`:
  
```r
# Install scDEDM from GitHub
devtools::install_github("hth-buaa/scDEDM2")

# Load the package
library(scDEDM)
```

## Quick Start

### Data Requirements
To use scDEDM, you need:
1. **Cell-matched scRNA-seq and scATAC-seq data**: Both datasets must have identical cell names.
2. **Cell type annotations**: A metadata column named `seurat_annotations` in both Seurat objects, containing matching cell type labels for each cell.
3. **scATAC-seq fragment file**: A `_fragments.tsv.gz` file and its corresponding `.tbi` index file (required for Signac-based analysis).

### Full Tutorial Link
For a complete step-by-step tutorial with code examples, please refer to the scDEDM2_guide_line.R or scDEDM2_guide_line.txt. This README focuses on the workflow logic and key parameters; the tutorial provides executable code for all steps.

## Detailed Workflow

### Part 0: Data Loading and Preprocessing
This initial step prepares your input data for downstream analysis by ensuring data quality, consistency, and compatibility with scDEDM.

#### Loading Seurat Objects
- **Example Dataset**: scDEDM provides a workflow using the `pbmcMultiome` dataset (10x Genomics PBMC multi-ome data) in R package SeuratData. 
- **Custom Data**: For your own data, load cell-matched scRNA-seq and scATAC-seq Seurat objects. If data comes from separate experiments, first establish one-to-one cell correspondence and unify cell identifiers.

#### Upgrading to Seurat v5
scDEDM is built on Seurat v5. Convert your scRNA-seq assay from Seurat v4 `Assay` to Seurat v5 `Assay5` if needed.

#### Filtering Abnormal Chromosomes
Remove peaks from non-standard chromosomes (e.g., contigs) to focus on biologically relevant regions.

#### Adding Cell Labels
Ensure both Seurat objects have a `seurat_annotations` column in their metadata with matching cell type labels. If missing, manually add or adjust this column.

#### Sampling Cells
- **Filter Cell Types**: Subset data to focus on cell types of interest (e.g., monocytes in the example).
- **Unify Labels**: For related cell types (e.g., CD14 Mono and CD16 Mono), you need to unify their labels (e.g., "Mono") to analyze their combined differentiation trajectory.

#### Downsampling Genes and Peaks (Optional)
To speed up testing, you can downsample data while preserving biological signal:
- **Cells**: Retain top cells by combined non-zero frequency of RNA expression and chromatin accessibility.
- **Genes**: Filter out non-coding RNAs (e.g., AC/AL/LINC prefixes) and mitochondrial genes, then select top variable genes.
- **Peaks**: Select top variable peaks based on accessibility variability.

#### Preparing Fragment Files
scATAC-seq analysis requires a `_fragments.tsv.gz` file and `.tbi` index:
- **Download**: If provided, download the fragment file from the source (e.g., 10x Genomics).
- **Generate**: If missing, generate the fragment file from the scATAC-seq count matrix and compress it with `bgzip`, then create an index with `tabix` (Linux tools).

### Part 1: Data Preprocessing
This step extracts key regulatory elements (TGs and TFs) and generates matrices for GRN construction.

#### Annotating Peaks and Identifying Target Genes (TGs)
Use `ChIPseeker` to annotate chromatin accessibility peaks and identify TGs with open promoters:
- **Promoter Range**: Define a promoter search range (e.g., ±3000 bp around TSS; smaller ranges reduce runtime but may miss distal elements).
- **Output**: Annotation plots, peak annotation tables, and a list of TGs.

#### Identifying Transcription Factors (TFs)
Extract TFs present in your scRNA-seq data from the JASPAR2024 database:
- **Species**: Set to your organism (e.g., "Homo sapiens").
- **Collections**: Include JASPAR collections (e.g., "CORE", "PHYLOFACTS") to expand TF coverage.

#### Generating Expression and Activity Matrices
Process multi-omics data to extract:
- **TF/TG Expression Matrices**: From scRNA-seq data.
- **Gene Activity Matrices**: From scATAC-seq data (using Signac).
- **Input**: Provide the path to your fragment file (with `.tbi` index).

### Part 2: Pseudotime Inference and Cell Grouping
This step models cell differentiation trajectories and groups cells to handle data sparsity.

#### Pseudotime Analysis and Branch Partitioning
- **Quality Control**: Filter genes and cells for the selected cell type(s) to ensure data quality.
- **Trajectory Inference**: Use Monocle2 (with DDRTree) to infer pseudotemporal ordering of cells.
- **Branch Partitioning**: Interactively define branches based on trajectory visualization (e.g., 2 branches for monocyte differentiation).

#### Cell Grouping Along Pseudotime
- **Extract Temporal Profiles**: For each gene, extract expression (TFs/TGs) and activity (TGs) values along pseudotime.
- **Filter Genes/Cells**: Use `alpha_gene` (max missing rate for genes, e.g., 0.7) and `alpha_cell` (max missing rate for cells, e.g., 0.9) to remove low-quality data.
- **Adaptive Cell Grouping**: Group cells along pseudotime using a logarithmic fitting model (where xrepresents the number of cells in a branch, and yrepresents the minimum non-zero count of gene expression and activity values across all genes within each cell group for that branch) to ensure each group has sufficient non-zero values. Define fitting points (e.g., (200,5), (500,10), (1000,15)) to balance group size and number.

### Part 3: Constructing Initial Gene Regulatory Networks (iGRNs)
This step builds a preliminary GRN using TFBS predictions.

#### TFBS PWM-Based Initial GRN Construction
- **PWM Matching**: Use JASPAR2024 PWMs to predict TF binding sites in promoter regions of TGs.
- **Min Score**: Set `min_score_for_matchPWM` (e.g., "0%") to retain all potential binding sites (stringent filtering occurs later).
- **Output**: A TF-TG association matrix with initial regulatory strength scores (theta_i).

#### Setting Regulatory Threshold (Tao)
Calculate a cell-type-specific threshold `tao` to binarize the iGRN:
- **Definition**: Regulatory relationships with theta_i > tao are considered significant.
- **Calculation**: Tao is set to the minimum non-zero theta_i minus a small constant (lower bound: 0.005).

#### Branch-Specific iGRNs
Extract subnetworks for each branch by filtering the global iGRN to include only TFs and TGs present in that branch.

### Part 4: Training DEDM-Based Predictive Models
This step trains regression models to predict regulatory strengths using a hybrid optimization approach.

#### Training Set Partitioning
Select high-confidence training samples from the iGRN:
- **Top Theta Selection**: For each TF and TG, retain top interactions by theta_i.
- **Filter Thresholds**: Apply quantile and absolute thresholds to refine the training set.
- **Refinement (Optional)**: Require mutual top membership (TF in TG's top regulons and vice versa) to further refine training data.

#### Model Training with Hybrid Optimization
Train models for each branch using a combination of genetic algorithms (to escape local optima) and gradient ascent (to speed convergence):
- **Key Parameters**:
  - `max_epochs`: Maximum training rounds (e.g., 1000; training stops early if convergence is reached).
  - `early_stop`: Stop if no improvement in N epochs (e.g., 5; balances performance and overfitting).
  - `popSize_ga`/`maxiter_ga`: Genetic algorithm population size (e.g., 50) and generations (e.g., 150).
  - `iterations_grad`: Gradient ascent steps per iteration (e.g., 30).
  - `theta_difference_threshold`/`fit_loss_threshold`: Quality thresholds for accepting models (e.g., 0.1 each).
- **Parallel Computing**: Use multiple cores (set `ncores` to available cores minus 1; Linux recommended).

### Part 5: Predicting and Refining GRNs
This step uses trained models to infer high-confidence GRNs (TF-TG) and link them to chromatin peaks.

#### Preparing Data for GRN Inference
Generate all possible TF-TG pairs and attach their branch-specific cell group state features (TF expression, TG activity, TG expression).

#### Inferring Branch-Specific GRNs
- **Predict Regulatory Strengths**: Apply trained models to predict regulatory strength (theta_p) and fitting loss for all TF-TG pairs.
- **Parallel Processing**: Partition data to distribute load across cores (set `n_part` and `each` to balance workload).

#### Selecting Best Predictions
Filter predictions to retain high-quality interactions:
- **Initial Thresholds**: Use quantile (e.g., 0.8 for theta_p, 0.2 for loss) and absolute bounds (e.g., 0.2 for theta_p, 0.2 for loss).
- **Branch Integration**: Aggregate results across branches using weighted averaging (by cell count) to get a unified GRN (TF-TG).

#### Filtering GRNs (Loose and Strict)
Apply a second round of filtering to generate two GRN versions:
- **Loose GRN**: Retains regulatory relationships present in any branch (comprehensive).
- **Strict GRN**: Retains only relationships present in all branches (core).
- **Recommendation**: Use the loose GRN for most analyses (retains more biological signal).

#### Inferring Peak-TG and TF-Peak Relationships
Link the GRN to chromatin accessibility peaks:
- **Peak-TG Inference**:
  - Identify peaks in promoter regions of TGs in the filtered GRN.
  - Calculate importance scores for peak-TG pairs using features like accessibility, expression, and non-zero ratios.
  - Filter pairs using a weighted composite score (customize weights based on your research focus).
- **TF-Peak Inference**:
  - Link TFs in the GRN to peaks associated with their target genes.
  - Calculate importance scores for TF-peak pairs (include theta_p, theta_i, and fit_loss as features).
  - Filter pairs, ensuring at least one peak per TF is retained.

#### Constructing TF-Peak-TG eGRNs
Integrate TF-peak and peak-TG predictions to build the final eGRN (TF-peak-TG) :
- **Regulation Types**: Classify relationships as "TF-peak-TG" (both components supported), "TF-peak" (only TF-peak supported), or "peak-TG" (only peak-TG supported).
- **Output**: A full eGRN with all features and a simplified version with key columns (TF, peak, TG, regulation type).

## Output Description
scDEDM generates the following key outputs:
1. **Processed Data**: Filtered Seurat objects, expression/activity matrices, and pseudotime information.
2. **Initial GRNs (iGRNs)**: Global and branch-specific TF-TG networks with initial regulatory scores.
3. **Trained Models**: Branch-specific regression models for regulatory strength prediction.
4. **Filtered Predictive GRNs (pGRN)**: Loose and strict versions of the TF-TG GRN with predicted regulatory strengths.
5. **Peak-TG/TF-Peak Predictions**: Filtered relationships with importance scores.
6. **Final Predictive eGRN (peGRN)**: Integrated TF-peak-TG network with regulation type annotations.

All intermediate results can be saved using `saveRDS` (or `qs::qsave` for large files) to enable step-by-step reanalysis.

## Memory Management Tips
- **Parallel Computing**: scDEDM uses parallel processing. If memory is exhausted, reduce the number of cores (`ncores`).
- **Restart R**: After parallel steps, memory may not be fully released. Restart R and reload intermediate results to free memory.
- **Independent Functions**: All scDEDM functions can run independently. Save intermediate results after each step to avoid re-running the entire workflow.

## Common Issues and Solutions
1. **Data Download Failure**:
   - For `pbmcMultiome`, retry downloading when network speed is good.
   - If persistent, manually download the data tarball and install it locally.

2. **Memory Errors**:
   - Reduce `ncores` or downsample data (cells/genes/peaks) for testing.
   - Use `qs::qsave`/`qs::qread` instead of `saveRDS`/`readRDS` for faster I/O and smaller file sizes.

3. **Unsatisfactory Cell Grouping**:
   - Adjust `points_x_for_fitting_nc_and_nmin` and `points_y_for_fitting_nc_and_nmin` to change the relationship between total cells and minimum group size.
   - Check the number of groups and minimum non-zero counts per group to ensure statistical power.

4. **Model Training Convergence Issues**:
   - Increase `max_epochs` or `early_stop` to allow longer training.
   - Adjust `theta_difference_threshold` or `fit_loss_threshold` to relax quality requirements.
   - Increase `popSize_ga` or `maxiter_ga` to improve the genetic algorithm's ability to escape local optima.
                                                            
## Citation
If you use scDEDM in your research, please cite:
[...]
                                                            
## Related Links
1. R package is in https://github.com/hth-buaa/scDEDM2.
2. Specific Guidelines is in https://github.com/hth-buaa/scDEDM2/blob/main/scDEDM2_guide_line.R or https://github.com/hth-buaa/scDEDM2/blob/main/scDEDM2_guide_line.txt.
3. Experiment results of the paper is in ... . 
4. The code for figures and some mediate analysis results are in ... .
