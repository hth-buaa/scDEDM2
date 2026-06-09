#' @title Extract Gene Expression and Activity Profiles Along Pseudotime
#'
#' @description
#' This function processes single-cell multi-omics data to extract temporal expression and activity profiles for transcription factors (TFs) and target genes (TGs) along developmental trajectories.
#' It performs rigorous quality filtering, maps cells to pseudotime indices, normalizes data, and returns organized temporal profiles for downstream GRN inference.
#' This function processes single-cell multi-omics data to extract temporal expression and activity profiles for transcription factors (TFs) and target genes (TGs) along developmental trajectories.
#' It performs rigorous quality filtering based on multi-dimensional gene importance scores, maps cells to pseudotime indices, normalizes data, and returns organized temporal profiles for downstream GRN inference.
#'
#' @param interest_cell_type_Branches The output of function order_pseudotime_and_divide_branches.
#' @param interest_cell_type_data The output of function get_interest_cell_type_data
#' @param alpha_TF Numeric (0 to 1). Threshold for TF filtering based on TF Importance Score (TFIS). TFs with TFIS above this quantile will be retained. Default is 0.7.
#' @param alpha_TG Numeric (0 to 1). Threshold for TG filtering based on TG Importance Score (TGIS). TGs with TGIS above this quantile will be retained. Default is 0.7.
#' @param alpha_cell Numeric (0 to 1). Threshold for cell filtering based on zero-expression rate. Cells with zero-expression rate higher than this value will be removed. Default is 0.9.
#' @param w_TFNER Numeric (-0.2 to 1). Weight for TF Non-Zero Expression Ratio in TFIS calculation. Default is 1.
#' @param w_TFNERR Numeric (-0.2 to 1). Weight for TF Non-zero Expression Ratio normalized Rank (max corresponds to 1, min to 0) in TFIS calculation. Default is 0.6.
#' @param w_TFNEA Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Expression Average in TFIS calculation. Default is 0.8.
#' @param w_TFNEAR Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Expression Average Rank (max corresponds to 1, min to 0) in TFIS calculation. Default is 0.4.
#' @param w_TFNAR Numeric (-0.2 to 1). Weight for TF Non-Zero Activity Ratio in TFIS calculation. Default is 0.9.
#' @param w_TFNARR Numeric (-0.2 to 1). Weight for TF Non-zero Activity Ratio normalized Rank (max corresponds to 1, min to 0) in TFIS calculation. Default is 0.5.
#' @param w_TFNAA Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Activity Average in TFIS calculation. Default is 0.7.
#' @param w_TFNAAR Numeric (-0.2 to 1). Weight for TF normalized (divided by max) Non-zero Activity Average Rank (max corresponds to 1, min to 0) in TFIS calculation. Default is 0.3.
#' @param w_TGNER Numeric (-0.2 to 1). Weight for TG Non-Zero Expression Ratio in TGIS calculation. Default is 0.3.
#' @param w_TGNERR Numeric (-0.2 to 1). Weight for TG Non-zero Expression Ratio normalized Rank (max corresponds to 1, min to 0) in TGIS calculation. Default is 0.1.
#' @param w_TGNEA Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Expression Average in TGIS calculation. Default is 0.3.
#' @param w_TGNEAR Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Expression Average Rank (max corresponds to 1, min to 0) in TGIS calculation. Default is 0.1.
#' @param w_TGNAR Numeric (-0.2 to 1). Weight for TG Non-Zero Activity Ratio in TGIS calculation. Default is 1.
#' @param w_TGNARR Numeric (-0.2 to 1). Weight for TG Non-zero Activity Ratio normalized Rank (max corresponds to 1, min to 0) in TGIS calculation. Default is 1.
#' @param w_TGNAA Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Activity Average in TGIS calculation. Default is 0.3.
#' @param w_TGNAAR Numeric (-0.2 to 1). Weight for TG normalized (divided by max) Non-zero Activity Average Rank (max corresponds to 1, min to 0) in TGIS calculation. Default is 0.1.
#' @param ncores See ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list structure organized by cell type, where for each cell type contains:
#' \itemize{
#' \item \code{Branches_TFE}: List of data frames containing temporal expression vectors for TFs in each branch
#' \item \code{Branches_TGA}: List of data frames containing temporal activity vectors for TGs in each branch
#' \item \code{Branches_TGE}: List of data frames containing temporal expression vectors for TGs in each branch
#' \item \code{Branches_n_cell}: List of cell counts for each branch
#' \item \code{Branches}: Original branch information with pseudotime mapping
#' \item \code{Branches_TFs}: List of TF names present in each branch
#' \item \code{Branches_TGs}: List of TG names present in each branch
#' }
#' All vectors are ordered by pseudotime indices within each branch.
#'
#' @details
#' The function performs the following processing steps for each cell type and branch:
#' \enumerate{
#' \item Calculates 8-dimensional features for both TFs and TGs, including non-zero expression ratio, its normalized rank, normalized non-zero expression average, its normalized rank, and the corresponding four metrics for activity data
#' \item The range for all weight parameters (those starting with w_) is recommended to be –0.2 to 1. Since all features are positive (higher feature values indicate greater importance), positive weights are expected. However, slightly negative weights (between –0.2 and 0) are acceptable for less important or irrelevant features.
#' \item Computes TF Importance Score (TFIS) and TG Importance Score (TGIS) as weighted sums of the 8 features using user-specified weights
#' \item Applies iterative quality filtering: retains TFs with TFIS > \code{alpha_TF} quantile, TGs with TGIS > \code{alpha_TG} quantile, and cells with zero-expression rate ≤ \code{alpha_cell}
#' \item Maps cells to pseudotime indices (k~ values) for temporal ordering
#' \item Normalizes both expression and activity matrices using Seurat's \code{NormalizeData} method
#' \item Extracts and orders temporal profiles for TFs and TGs
#' \item Returns organized data structures ready for dynamic GRN inference
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Alpha parameters are hyperparameters that significantly affect results
#' \item The 8-dimensional weighting scheme allows flexible integration of expression and activity data for gene selection
#' \item Output vectors are ordered by pseudotime, enabling time-series analysis
#' \item Uses parallel processing in Linux for efficient handling of multiple cell types
#' \item After the function runs, it is crucial to check the non-zero ratios in the resulting vectors: for all TF expression vectors, all TG activity vectors, and all TG expression vectors; excessively low ratios may adversely affect or prevent subsequent cell grouping
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_Branches = base::readRDS(paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_Branches.rds"
#' ))
#' interest_cell_type_data = base::readRDS(paste0(
#'   "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment",
#'   "/interest_cell_type_data.rds"
#' ))
#' alpha_TF = 0.7; alpha_TG = 0.7; alpha_cell = 0.9
#' w_TFNER = 1; w_TFNERR = 0.6; w_TFNEA = 0.8; w_TFNEAR = 0.4
#' w_TFNAR = 0.9; w_TFNARR = 0.5; w_TFNAA = 0.7; w_TFNAAR = 0.3
#' w_TGNER = 0.3; w_TGNERR = 0.1; w_TGNEA = 0.3; w_TGNEAR = 0.1
#' w_TGNAR = 1; w_TGNARR = 1; w_TGNAA = 0.3; w_TGNAAR = 0.1
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_TF = alpha_TF, alpha_TG = alpha_TG, alpha_cell = alpha_cell,
#'   w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
#'   w_TFNAR = w_TFNAR, w_TFNARR = w_TFNARR, w_TFNAA = w_TFNAA, w_TFNAAR = w_TFNAAR,
#'   w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
#'   w_TGNAR = w_TGNAR, w_TGNARR = w_TGNARR, w_TGNAA = w_TGNAA, w_TGNAAR = w_TGNAAR,
#'   ncores = ncores
#' )
#' # Check minimum nonzero frequency
#' min_nonzero_freq = 1.0
#' for (i in seq_along(interest_cell_type_genes_pseudotime_info)) {
#'   for (j in 1:3) {
#'     for (k in seq_along(interest_cell_type_genes_pseudotime_info[[i]][[j]])) {
#'       df = interest_cell_type_genes_pseudotime_info[[i]][[j]][[k]]
#'       nonzero_freq = rowSums(df != 0) / ncol(df)
#'       min_nonzero_freq = min(min_nonzero_freq, nonzero_freq)
#'     }
#'   }
#' }
#' cat("min_nonzero_freq: ", min_nonzero_freq, "\n"); rm(min_nonzero_freq, df, nonzero_freq)
#' base::saveRDS(
#'   interest_cell_type_genes_pseudotime_info,
#'   file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds"
#' )
#' }
get_genes_pseudotime_info = function(
    interest_cell_type_Branches = interest_cell_type_Branches,
    interest_cell_type_data = interest_cell_type_data,
    alpha_TF = 0.7, alpha_TG = 0.7, alpha_cell = 0.9,
    w_TFNER = 1, w_TFNERR = 0.6, w_TFNEA = 0.8, w_TFNEAR = 0.4,
    w_TFNAR = 0.9, w_TFNARR = 0.5, w_TFNAA = 0.7, w_TFNAAR = 0.3,
    w_TGNER = 0.3, w_TGNERR = 0.1, w_TGNEA = 0.3, w_TGNEAR = 0.1,
    w_TGNAR = 1, w_TGNARR = 1, w_TGNAA = 0.3, w_TGNAAR = 0.1,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Obtaining Gene Expression/Activity Vectors along Pseudotime ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.2 Data Processing - Cell Grouping"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining a function to process each cell type individually.
  get_genes_pseudotime_info_process = function(
    Branches,
    interest_scRNAseq_for_GRN, interest_scATACseq_for_GRN,
    interest_gene_for_GRN, interest_TGs, interest_TFs,
    alpha_TF, w_TFNER, w_TFNERR, w_TFNEA, w_TFNEAR, w_TFNAR, w_TFNARR, w_TFNAA, w_TFNAAR,
    alpha_TG, w_TGNER, w_TGNERR, w_TGNEA, w_TGNEAR, w_TGNAR, w_TGNARR, w_TGNAA, w_TGNAAR,
    alpha_cell)
  {
    # Obtaining data for each branch.
    message("Initializing lists for storing scRNAseq, scATACseq, gene_for_GRN, TGs, TFs, n_cell of each branch.")
    Branches_RNA = base::list()
    Branches_ATAC = base::list()
    Branches_gene_for_GRN = base::list()
    Branches_TGs = base::list()
    Branches_TFs = base::list()
    Branches_n_cell = base::list()

    # Define feature calculation function
    calculate_features = function(exp_mat, # Gene expression matrix (rows: genes, columns: cells)
                                  act_mat, # Gene activity matrix (rows: genes, columns: cells)
                                  genes, # Vector of gene names (e.g., TF or TG names)
                                  weights, # Numeric vector of weights for the 8 features (order: NER, NERR, NEA, NEAR, NAR, NARR, NAA, NAAR)
                                  prefix) { # Prefix string for column names of the output data frame (e.g., "TF" or "TG")
      # Calculate non-zero rate (NER/NAR)
      ner = base::rowMeans(exp_mat != 0)
      nar = base::rowMeans(act_mat != 0)

      # Calculate average of non-zero values (NEA/NAA) (replace 0 with NA then use rowMeans)
      exp_na = exp_mat
      exp_na[exp_na == 0] = NA
      nea = base::rowMeans(exp_na, na.rm = TRUE)
      rm(exp_na)

      act_na = act_mat
      act_na[act_na == 0] = NA
      naa = base::rowMeans(act_na, na.rm = TRUE)
      rm(act_na)

      # Handle NA by replacing with 0
      ner[base::is.na(ner)] = 0
      nea[base::is.na(nea)] = 0
      nar[base::is.na(nar)] = 0
      naa[base::is.na(naa)] = 0

      # Normalize (divide by maximum value)
      ner = ner / max(ner)
      nea = nea / max(nea)
      nar = nar / max(nar)
      naa = naa / max(naa)

      # Calculate ranks
      n = length(genes)
      nerr = base::rank(ner) / n
      near = base::rank(nea) / n
      narr = base::rank(nar) / n
      naar = base::rank(naa) / n

      # Calculate importance score
      feat_mat = base::cbind(ner, nerr, nea, near, nar, narr, naa, naar)
      is_score = base::as.vector(feat_mat %*% weights)

      # Output feature matrix
      df = data.frame(
        feat_mat,
        IS = is_score,
        row.names = genes,
        check.names = FALSE
      )
      base::colnames(df) = paste0(prefix, c("NER", "NERR", "NEA", "NEAR", "NAR", "NARR", "NAA", "NAAR", "IS"))

      return(df)
    }

    for (i in base::seq_along(Branches)) {
      message("Obtaining gene expression count matrix for branch ", i, ".")
      Branches_RNA[[i]] = interest_scRNAseq_for_GRN[, base::rownames(Branches[[i]])]
      rna_assay = SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")

      message("Extracting gene activity count matrix for branch ", i, ".")
      Branches_ATAC[[i]] = interest_scATACseq_for_GRN[base::rownames(Branches_RNA[[i]]), base::colnames(Branches_RNA[[i]])]
      atac_assay = SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "counts")

      message("Calculating feature matrix for TFs in branch ", i, ".")
      TF_features = calculate_features(
        exp_mat = base::as.matrix(rna_assay[interest_TFs, ]),
        act_mat = base::as.matrix(atac_assay[interest_TFs, ]),
        genes = interest_TFs,
        weights = c(w_TFNER, w_TFNERR, w_TFNEA, w_TFNEAR, w_TFNAR, w_TFNARR, w_TFNAA, w_TFNAAR),
        prefix = "TF"
      )

      message("Calculating feature matrix for TGs in branch ", i, ".")
      TG_features = calculate_features(
        exp_mat = base::as.matrix(rna_assay[interest_TGs, ]),
        act_mat = base::as.matrix(atac_assay[interest_TGs, ]),
        genes = interest_TGs,
        weights = c(w_TGNER, w_TGNERR, w_TGNEA, w_TGNEAR, w_TGNAR, w_TGNARR, w_TGNAA, w_TGNAAR),
        prefix = "TG"
      )
      rm(rna_assay, atac_assay)

      message("Gene filtering for branch ", i, ".")
      TFs_filter = base::rownames(TF_features)[TF_features$TFIS > stats::quantile(TF_features$TFIS, alpha_TF)]
      TGs_filter = base::rownames(TG_features)[TG_features$TGIS > stats::quantile(TG_features$TGIS, alpha_TG)]
      genes_filter = base::sort(base::unique(c(TFs_filter, TGs_filter)))
      Branches_RNA[[i]] = subset(Branches_RNA[[i]], features = genes_filter)
      Branches_ATAC[[i]] = subset(Branches_ATAC[[i]], features = genes_filter)
      rm(genes_filter)

      message("Cell filtering for branch ", i, ".")
      zero_freq_RNA = base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")) == 0))
      zero_freq_ATAC = base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "counts")) == 0))
      cells_filter = base::intersect(
        base::names(zero_freq_RNA[zero_freq_RNA <= alpha_cell]),
        base::names(zero_freq_ATAC[zero_freq_ATAC <= alpha_cell])
      )
      rm(zero_freq_RNA, zero_freq_ATAC)
      Branches_RNA[[i]] = subset(Branches_RNA[[i]], cells = cells_filter)
      Branches_ATAC[[i]] = subset(Branches_ATAC[[i]], cells = cells_filter)
      rm(cells_filter)

      message("Retrieving TGs and TFs for branch ", i, ".")
      Branches_gene_for_GRN[[i]] = base::intersect(interest_gene_for_GRN, base::rownames(Branches_RNA[[i]]))
      Branches_TGs[[i]] = base::intersect(base::intersect(interest_TGs, Branches_gene_for_GRN[[i]]), TGs_filter); rm(TGs_filter)
      Branches_TFs[[i]] = base::intersect(base::intersect(interest_TFs, Branches_gene_for_GRN[[i]]), TFs_filter); rm(TFs_filter)

      message("Retrieving the cell number for branch ", i, ".")
      Branches_gene_for_GRN[[i]] = base::union(Branches_TGs[[i]], Branches_TFs[[i]])
      Branches_RNA[[i]] = Branches_RNA[[i]][base::intersect(base::rownames(Branches_RNA[[i]]), interest_gene_for_GRN), ]
      Branches_ATAC[[i]] = Branches_ATAC[[i]][base::rownames(Branches_RNA[[i]]), ]
      Branches_n_cell[[i]] = base::ncol(Branches_RNA[[i]])
    }
    rm(calculate_features)

    # Mapping cells in expression and activity matrices to pseudotime point indices k~ for each branch.
    message("Mapping cells in expression and activity matrices to pseudotime point indices k~ for each branch.")
    for (i in base::seq_along(Branches)) {
      message("Creating a cell-to-pseudotime mapping for branch ", i, ".")
      Branches[[i]]["cell"] = base::rownames(Branches[[i]])
      Branches[[i]]["k~"] = 1:base::nrow(Branches[[i]])
      name_mapping = stats::setNames(Branches[[i]]$`k~`, Branches[[i]]$cell)

      message("Renaming cells in branch ", i, " to pseudotime indices for each branch.")
      base::colnames(Branches_RNA[[i]]) = name_mapping[base::colnames(Branches_RNA[[i]])]
      base::colnames(Branches_ATAC[[i]]) = name_mapping[base::colnames(Branches_ATAC[[i]])]
    }

    # Normalizing expression and activity matrices for each branch.
    message("Normalizing expression and activity matrices for each branch.")
    for (i in base::seq_along(Branches)) {
      message("Normalizing scRNA-seq matrix for branch ", i, ".")
      Branches_RNA[[i]] = Seurat::NormalizeData(Branches_RNA[[i]])

      message("Normalizing scATAC-seq matrix for branch ", i, ".")
      Branches_ATAC[[i]] = Seurat::NormalizeData(Branches_ATAC[[i]])
    }

    # For each branch, obtaining temporal expression vectors for each TF, temporal expression vectors for each TG, and temporal activity vectors for each TG.
    message("Initializing lists for storing temporal expression vectors of TFs, temporal expression vectors of TGs, temporal activity vectors of TGs across all branches.")
    Branches_TFE = base::list()
    Branches_TGA = base::list()
    Branches_TGE = base::list()
    for (i in base::seq_along(Branches)) {
      message("Obtaining temporal expression vectors for each TF in branch ", i, ".")
      Branches_TFE[[i]] =
        base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data")))[Branches_TFs[[i]], ]
      Branches_TFE[[i]] = Branches_TFE[[i]][, base::order(base::as.integer(base::colnames(Branches_TFE[[i]])))]

      message("Obtaining temporal activity vectors for each TG in branch ", i, ".")
      Branches_TGA[[i]] =
        base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "data")))[Branches_TGs[[i]], ]
      Branches_TGA[[i]] = Branches_TGA[[i]][, base::order(base::as.integer(base::colnames(Branches_TGA[[i]])))]

      message("Obtaining temporal expression vectors for each TG in branch ", i, ".")
      Branches_TGE[[i]] =
        base::as.data.frame(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data"))[Branches_TGs[[i]], ]
      Branches_TGE[[i]] = Branches_TGE[[i]][, base::order(base::as.integer(base::colnames(Branches_TGE[[i]])))]
    }

    return(base::list("Branches_TFE" = Branches_TFE,
                      "Branches_TGA" = Branches_TGA,
                      "Branches_TGE" = Branches_TGE,
                      "Branches_n_cell" = Branches_n_cell,
                      "Branches" = Branches,
                      "Branches_TFs" = Branches_TFs,
                      "Branches_TGs" = Branches_TGs))
  }

  ### Invoking the newly defined function in parallel.
  interest_cell_type_genes_pseudotime_info = parallel::mclapply(
    X = base::names(interest_cell_type_Branches),
    FUN = function(cell_type) {
      message("Preparing to extract pseudotime vectors of TFE, TGA, and TGE across all branches in cell type ", cell_type, ".")
      get_genes_pseudotime_info_process(
        Branches = interest_cell_type_Branches[[cell_type]],
        interest_scRNAseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scRNAseq_for_GRN"]],
        interest_scATACseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scATACseq_for_GRN"]],
        interest_gene_for_GRN = interest_cell_type_data[[cell_type]][["interest_gene_for_GRN"]],
        interest_TGs = interest_cell_type_data[[cell_type]][["interest_TGs"]],
        interest_TFs = interest_cell_type_data[[cell_type]][["interest_TFs"]],
        alpha_TF = alpha_TF,
        w_TFNER = w_TFNER, w_TFNERR = w_TFNERR, w_TFNEA = w_TFNEA, w_TFNEAR = w_TFNEAR,
        w_TFNAR = w_TFNAR, w_TFNARR = w_TFNARR, w_TFNAA = w_TFNAA, w_TFNAAR = w_TFNAAR,
        alpha_TG = alpha_TG,
        w_TGNER = w_TGNER, w_TGNERR = w_TGNERR, w_TGNEA = w_TGNEA, w_TGNEAR = w_TGNEAR,
        w_TGNAR = w_TGNAR, w_TGNARR = w_TGNARR, w_TGNAA = w_TGNAA, w_TGNAAR = w_TGNAAR,
        alpha_cell = alpha_cell
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_genes_pseudotime_info) = base::names(interest_cell_type_Branches)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Obtaining Gene Expression/Activity Vectors along Pseudotime ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_genes_pseudotime_info)
}
