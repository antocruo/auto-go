#' @title  Differential Gene Expression Analysis
#'
#' @description This function allows to perform a differential gene expression analysis using the DESeq2 package.
#' @description The principal DESeq2 workflow is employed. Raw counts are rounded, if it is necessary, because integers are needed for DESeq2 to run. In case the user provides an .rds file, the tool makes sure that in the assay of the SummarizedExperiment only counts are stored. Then the function DESeqDataSetFromMatrix (or DESeqDataSet for .rds) is employed. A prefiltering is applied in order to remove all genes having sum along the subjects less than 10. The Differential Expression Analysis is performed by employing the function DE and the normalized data matrix is stored in the working directory. This will be useful for further analysis and visualizations. Results are extracted for each comparison of interest. Subfolders will be generated in the "outfolder" (default: ./results) with the name of the comparisons made. NOTE: as standard we use the nomenclature "CONTROL_vs_TREATMENT", i.e. the control is on the left. Inside each comparison subfolder we will find a .tsv file with the complete differential analysis and other subfolder based on the pvalue and the log2FC thresholds; inside this we will find a .tsv file with the results of the only filtered genes. This subfolders will be divided in other subfolders as "up_genes", "down_genes" and "up_down_genes". Look the path flow chart at the end of this tutorial (Figure 1).
#' @param counts The path to raw counts file. Accepted file formats are tab or comma-separated files (.tsv, .csv), .txt files, .rds. Genes must be on rows, samples on columns.
#' @param groups Sample information table needed by DESeq2 (e.g. 'colData'). A data frame with at least two columns: one for samples, one for a grouping variable (See examples).
#' @param comparisons Table of comparisons based on the grouping variable in 'groups' table (See examples). It should be a data.frame with column 'treatment' and column 'control'. It is possible to provide the path to a .txt file.
#' @param padj_threshold Threshold value for adjusted p-value significance (Defaults to 0.05).
#' @param log2FC_threshold Threshold value for log2(Fold Change) for considering genes as differentially expressed.
#' @param pre_filtering Removes genes which sum in the raw counts is less than 10 (Default = TRUE).
#' @param save_excel Allows to save all the output tables in .xlsx format (Default = FALSE).
#' @param outfolder The name to assign to the folder for output saving. (Default = "./results").
#' @param del_csv Specify the delimiter of the .csv file, default is ",". This is because opening .csv files with Excel messes up the format and changes the delimiter in ";".
#' @param calculate_tpm Logical. If TRUE, TPM values will be calculated and added to results (Default = FALSE).
#' @param gene_length_file Path to gene length file. If NULL and calculate_tpm = TRUE, the default gene_length.txt file from the package will be used. Tab-delimited file with gene names in first column and lengths in second column.
#' @return No return value. Files will be produced as part of normal execution.
#' @examples
#' sample <- c("Pat_1", "Pat_2", "Pat_3", "Pat_4", "Pat_5", "Pat_6")
#' group <- c("CTRL", "CTRL", "TREAT_A", "TREAT_A", "TREAT_B", "TREAT_B")
#' groups <- data.frame(sample, group)
#' treatment <- c("TREAT_A", "TREAT_B", "TREAT_A")
#' control <- c("CTRL", "CTRL", "TREAT_B")
#' comparisons <- data.frame(treatment, control)
#' \dontrun{
#' deseq_analysis(counts,
#'   groups,
#'   comparisons,
#'   padj_threshold = 0.05,
#'   log2FC_threshold = 0,
#'   pre_filtering = T,
#'   save_excel = F,
#'   outfolder = "./results",
#'   del_csv = ",",
#'   calculate_tpm = TRUE,
#'   gene_length_file = "gene_length.txt",
#' )
#' }
#' @export


deseq_analysis <- function(counts,
                           groups,
                           comparisons,
                           padj_threshold = 0.05,
                           log2FC_threshold = 0,
                           pre_filtering = TRUE,
                           save_excel = FALSE,
                           outfolder = "./results",
                           del_csv = ",",
                           calculate_tpm = FALSE,
                           gene_length_file = NULL) {
  
  # function to calculate TPM
  calculate_tpm_values <- function(counts_matrix, gene_lengths) {
    
    # match gene lengths to counts matrix
    matched_lengths <- gene_lengths[match(rownames(counts_matrix), gene_lengths[[1]]), 2]
    
    n_total <- length(matched_lengths)
    n_missing <- sum(is.na(matched_lengths))
    n_found <- n_total - n_missing
    
    if (n_missing > 0) {
      pct_missing <- round(n_missing / n_total * 100, 1)
      pct_found <- round(n_found / n_total * 100, 1)
      message(sprintf("TPM matching: %d/%d genes (%.1f%%) matched successfully", 
                      n_found, n_total, pct_found))
      message(sprintf("TPM matching: %d/%d genes (%.1f%%) not found in length file", 
                      n_missing, n_total, pct_missing))
      
      # show some missing genes
      missing_genes <- rownames(counts_matrix)[is.na(matched_lengths)]
      message(sprintf("Example missing genes: %s", paste(head(missing_genes, 5), collapse = ", ")))
      
      valid_idx <- !is.na(matched_lengths)
      counts_matrix <- counts_matrix[valid_idx, ]
      matched_lengths <- matched_lengths[valid_idx]
    } else {
      message(sprintf("TPM matching: All %d genes matched successfully!", n_total))
    }
    
    # Calculate RPK 
    rpk <- counts_matrix / (matched_lengths / 1000)
    
    # Calculate TPM
    tpm <- sweep(rpk, 2, colSums(rpk) / 1e6, FUN = "/")
    
    return(tpm)
  }
  
  # Check TPM requirements
  if (calculate_tpm) {
    if (is.null(gene_length_file)) {
      # gene length file from package
      gene_length_file <- system.file("extdata", "gene_length.txt", package = "autoGO")
      if (gene_length_file == "") {
        stop("gene_length_file must be provided when calculate_tpm = TRUE, or a default 'gene_length.txt' file must be available in the package.")
      }
      message("Using default gene length file from package: gene_length.txt")
    }
    if (!file.exists(gene_length_file)) {
      stop("Gene length file not found: ", gene_length_file)
    }
    gene_lengths <- read_delim(gene_length_file, delim = "\t", col_types = cols())
    if (ncol(gene_lengths) < 2) {
      stop("Gene length file must have at least 2 columns (e.g., gene_name and length)")
    }
    # Rename columns to standard names for consistency
    colnames(gene_lengths)[1:2] <- c("gene_id", "length")
  
    message(sprintf("Gene length file loaded: %d genes", nrow(gene_lengths)))
  }
  
  if (grepl(".tsv", counts)[1]) {
    counts <- read_delim(counts, col_types = cols(), delim = "\t")
  } else if (grepl(".csv", counts)[1]) {
    counts <- read_delim(counts, col_types = cols(), delim = del_csv)
  } else if (grepl(".rds", counts)[1]) {
    counts <- readRDS(counts)
  } else if (grepl(".txt", counts)[1]) {
    counts <- read_delim(counts, delim = "\t", col_types = cols())
  } else if (is.data.frame(counts)) {
    counts <- counts
  } else {
    warning("Provide a file .tsv, .csv, .rds or .txt tab-separated, or directly a data.frame")
  }

  if (!is.data.frame(groups)[1] & grepl(".txt", groups)[1]) {
    groups <- read_delim(groups, delim = "\t", col_types = cols())
  } else if (dim(groups)[2] < 2) {
    warning("Please provide a two column file as groups list as .txt file or data.frame.")
  }
  if (!is.data.frame(comparisons)[1] & grepl(".txt", comparisons)[1]) {
    comparisons <- read_delim(comparisons, delim = "\t", col_types = cols())
  } else if (dim(comparisons)[2] < 2) {
    warning("Please provide a two column file as comparison list as .txt file or data.frame.")
  }

  if (class(counts)[1] == "SummarizedExperiment") {
    SummarizedExperiment::assays(counts) <- SummarizedExperiment::assay(counts, 1)
    SummarizedExperiment::assay(counts) <- as.matrix(round(SummarizedExperiment::assay(counts)))

    rownames(groups) <- groups[[1]]
    if (!all(colnames(counts) == groups[[1]])) {
      groups <- groups[match(groups[[1]], colnames(counts)), ]
    }

    colData(counts) <- data.frame(groups)
    dds <- DESeqDataSet(se = counts, design = ~group)
  } else {
    colnames(counts)[1] <- "gene_id"
    counts <- counts %>%
      dplyr::distinct(.data$gene_id, .keep_all = TRUE) %>% # removing duplicated genes
      textshape::column_to_rownames(loc = "gene_id") %>%
      dplyr::relocate(sort(tidyselect::peek_vars())) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round))

    if (dim(counts)[2] != dim(groups)[1]) {
      stop("Please provide all counts columns in your groups dataframe.")
    }

    if (!all(colnames(counts) == groups[[1]])) {
      groups <- groups[match(colnames(counts), groups[[1]]), ]
    }
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = groups, design = ~group)
  }

  # Store raw counts before filtering for TPM calculation
  raw_counts <- counts(dds)
  
  # Diagnostic info for TPM
  if (calculate_tpm) {
    message(sprintf("Counts matrix has %d genes", nrow(raw_counts)))
  }
  
  # pre-filtering
  if (pre_filtering) {
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
  }

  # DE analysis
  dds <- DESeq(dds)
  cc <- counts(dds, normalized = T)
  vsd <- varianceStabilizingTransformation(dds, blind = T)
  vst <- SummarizedExperiment::assay(vsd)

  # Calculate TPM if requested
  tpm_matrix <- NULL
  if (calculate_tpm) {
    # Calculate TPM on raw counts (with optional gene mapping)
    tpm_matrix <- calculate_tpm_values(raw_counts, gene_lengths)
    
    # Filter TPM matrix to match filtered genes in dds
    # Keep only genes that are both in TPM matrix and in the filtered dds
    common_genes <- intersect(rownames(tpm_matrix), rownames(cc))
    tpm_matrix <- tpm_matrix[common_genes, ]
    
    n_missing_after_filter <- nrow(cc) - length(common_genes)
    if (n_missing_after_filter > 0) {
      pct_missing <- round(n_missing_after_filter / nrow(cc) * 100, 1)
      message(sprintf("After filtering: %d/%d genes (%.1f%%) in DESeq2 results lack TPM values (will show as NA).", 
                      n_missing_after_filter, nrow(cc), pct_missing))
    } else {
      message(sprintf("TPM successfully calculated for all %d filtered genes.", nrow(cc)))
    }
  }
  
  if (!dir.exists(outfolder)) {
    dir.create((outfolder), recursive = T)
  }
  write.table(cc, file.path(outfolder, "deseq_norm_data.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
  write.table(vst, file.path(outfolder, "deseq_vst_data.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
  
  if (calculate_tpm) {
    write.table(tpm_matrix, file.path(outfolder, "tpm_data.txt"), quote = F, sep = "\t", row.names = T, col.names = NA)
  }

  group_list <- split(groups, groups[[2]])
  for (i in seq(dim(comparisons)[1])) {
    c <- comparisons[i, ]

    a <- c[, 1] # treatment
    if (is.data.frame(a)) {
      a <- dplyr::pull(a)
    }
    b <- c[, 2] # control
    if (is.data.frame(b)) {
      b <- dplyr::pull(b)
    }

    a_group <- group_list[[a]][[1]]
    b_group <- group_list[[b]][[1]]

    rr <- results(dds, contrast = c("group", a, b), alpha = 0.05)
    rr <- na.omit(rr)
    rr <- as.data.frame(rr)

    means <- data.frame(mean_A = rowMeans(cc[, a_group]), mean_B = rowMeans(cc[, b_group]))
    names(means) <- c(paste0("mean_", a), paste0("mean_", b))
    means_vst <- data.frame(mean_A = rowMeans(vst[, a_group]), mean_B = rowMeans(vst[, b_group]))
    names(means_vst) <- c(paste0("mean_", a), paste0("mean_", b))

    rr <- merge(rr, means_vst, by = 0) %>%
      dplyr::rename(genes = .data$Row.names) %>%
      dplyr::select(-.data$baseMean, -.data$lfcSE, -stat) %>%
      dplyr::select(c(1:4, 6, 5)) %>%
      dplyr::arrange(.data$padj)

    # Add TPM means if calculated
    if (calculate_tpm && !is.null(tpm_matrix)) {
      # Get only genes present in TPM matrix
      genes_in_tpm <- intersect(rr$genes, rownames(tpm_matrix))
      
      if (length(genes_in_tpm) > 0) {
        means_tpm <- data.frame(
          gene = rownames(tpm_matrix),
          mean_tpm_A = rowMeans(tpm_matrix[, a_group, drop = FALSE]), 
          mean_tpm_B = rowMeans(tpm_matrix[, b_group, drop = FALSE])
        )
        names(means_tpm) <- c("genes", paste0("mean_tpm_", a), paste0("mean_tpm_", b))
        
        rr <- merge(rr, means_tpm, by = "genes", all.x = TRUE)
      }
    }
    
    filtered <- rr %>%
      dplyr::filter(rr$padj < padj_threshold & abs(rr$log2FoldChange) > log2FC_threshold)

    # generating folders
    groups_fold <- file.path(outfolder, paste0(b, "_vs_", a))
    groups_fold_filtered <- paste0("filtered_DE", "_thFC", log2FC_threshold, "_thPval", padj_threshold)
    groups_fold_filtered_path <- file.path(groups_fold, groups_fold_filtered)
    groups_fold_thresh_up_down <- file.path(groups_fold_filtered_path, "up_down_genes")
    groups_fold_thresh_up <- file.path(groups_fold_filtered_path, "up_genes")
    groups_fold_thresh_down <- file.path(groups_fold_filtered_path, "down_genes")

    # saving complete results
    if (!dir.exists(groups_fold_filtered_path)) dir.create(groups_fold_filtered_path, recursive = T)
    filename <- paste0("DE_", b, "_vs_", a, "_allres.tsv")
    write_tsv(rr, file.path(groups_fold, filename))

    if (save_excel) {
      filename <- paste0("DE_", b, "_vs_", a, "_allres.xlsx")
      openxlsx::write.xlsx(rr, file = file.path(groups_fold, filename), row.names = F)
    }

    # saving filtered results in different folders by thresholds
    filename <- paste0("filtered_DE_", b, "_vs_", a, "_thFC", log2FC_threshold, "_thPval", padj_threshold)
    write_tsv(filtered, file.path(groups_fold_filtered_path, paste0(filename, ".tsv")))
    if (save_excel) openxlsx::write.xlsx(filtered, file = file.path(groups_fold_filtered_path, paste0(filename, ".xlsx")), row.names = F)

    # saving gene lists
    if (!dir.exists(groups_fold_thresh_up_down)) dir.create(groups_fold_thresh_up_down, recursive = T)
    filename <- paste0("up_down_genes_list_", b, "_vs_", a, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".txt")
    write.table(filtered$genes, file.path(groups_fold_thresh_up_down, filename), quote = F, row.names = F, col.names = F)

    if (!dir.exists(groups_fold_thresh_up)) dir.create(groups_fold_thresh_up, recursive = T)
    filename <- paste0("up_genes_list_", b, "_vs_", a, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".txt")
    write.table(filtered$genes[filtered$log2FoldChange > 0], file.path(groups_fold_thresh_up, filename), quote = F, row.names = F, col.names = F)

    if (!dir.exists(groups_fold_thresh_down)) dir.create(groups_fold_thresh_down, recursive = T)
    filename <- paste0("down_genes_list_", b, "_vs_", a, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".txt")
    write.table(filtered$genes[filtered$log2FoldChange < 0], file.path(groups_fold_thresh_down, filename), quote = F, row.names = F, col.names = F)
  }
}