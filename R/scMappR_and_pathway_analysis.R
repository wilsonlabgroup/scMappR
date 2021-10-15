#' Generate cellWeighted_Foldchanges, visualize, and enrich.
#' 
#' This function generates cell weighted Fold-changes (cellWeighted_Foldchange), visualizes them in a heatmap, and completes pathway enrichment of cellWeighted_Foldchanges and the bulk gene list using g:ProfileR.
#'
#' This function generates cellWeighted_Foldchanges for every cell-type (see deconvolute_and_contextualize), as well as accompanying data such as cell-type proportions with the DeconRNA-seq, WGCNA, or DCQ methods.
#' Then, it generates heatmaps of all cellWeighted_Foldchanges, cellWeighted_Foldchanges overlapping with the signature matrix, the entire signature matrix, the cell-type preference values from the signature matrix that overlap with inputted differentially expressed genes.
#' Then, assuming there is available internet, it will complete gProfileR of the reordered cellWeighted_Foldchanges as well as a the ordered list of genes.
#' This function is a wrapper for deconvolute_and_contextualize and pathway_enrich_internal and the primary function within the package.
#' 
#' @rdname scMappR_and_pathway_analysis
#' @name scMappR_and_pathway_analysis
#' 
#' @param count_file Normalized (i.e. TPM, RPKM, CPM) RNA-seq count matrix where rows are gene symbols and columns are individuals. Inputted data should be a data.frame or matrix. A character vector to a tsv file where this data can be loaded is also acceptable. Gene symbols from the count file, signature matrix, and DEG list should all match (case sensitive, gene symbol or ensembl, etc.)
#' @param signature_matrix Signature matrix: a gene by cell-type matrix populated with the fold-change of gene expression in cell-type marker "i" vs all other cell-types. Object should be a data.frame or matrix.
#' @param DEG_list An object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted p-value, and the third the log2FC path to a .tsv file containing this info is also acceptable.
#' @param case_grep A character representing what designates the "cases" (i.e. upregulated is 'case' biased) in the columns of the count file. A numeric vector of the index of "cases" is also acceptable.  Tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases.
#' @param control_grep A character representing what designates the "control" (i.e. downregulated is 'control biased) in the columns of the count file. A numeric vector of the index of "control" is also acceptable.  Tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases.
#' @param max_proportion_change Maximum cell-type proportion change -- may be useful if there are many rare cell-type. Alternatively, if a cell-type is only present in one condition but not the other, it will prevent possible infinite or 0 cwFold-changes.
#' @param print_plots Whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed. The same name of the plots will be completed for top pathways.
#' @param plot_names The prefix of plot pdf files.
#' @param output_directory The name of the directory that will contain output of the analysis.
#' @param theSpecies human, mouse, or a species directly compatible with gProfileR (i.e. g:ProfileR). 
#' @param sig_matrix_size Maximum number of genes in signature matrix for cell-type deconvolution.
#' @param drop_unknown_celltype Whether or not to remove "unknown" cell-types from the signature matrix.
#' @param internet Whether you have stable Wifi (T/F).
#' @param up_and_downregulated Whether you are additionally splitting up/downregulated genes (T/F).
#' @param gene_label_size The size of the gene label on the plot.
#' @param number_genes The number of genes to cut-off for pathway analysis (good with many DEGs).
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' @param rda_path If downloaded, path to where data from scMappR_data is stored. 
#' @param newGprofiler Whether to use gProfileR or gprofiler2 (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#' @param deconMethod Which RNA-seq deconvolution method to use to estimate cell-type proporitons. Options are "WGCNA", "DCQ", or "DeconRNAseq"
#' 
#' 
#' @return List with the following elements:
#' \item{cellWeighted_Foldchanges}{Cellweighted Fold-changes for all differentially expressed genes.}
#' \item{paths}{Enriched biological pathways for each cell-type.}
#' \item{TFs}{Enriched TFs for each cell-type.}
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text geom_bar theme_classic xlab ylab scale_fill_manual element_line
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases ks.test dist shapiro.test mad
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#' @importFrom ADAPTS estCellPercent
#' @importFrom reshape melt
#'
#' @examples 
#' 
#' data(PBMC_example)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' toOut <- scMappR_and_pathway_analysis(count_file = bulk_normalized, 
#'                                       signature_matrix = odds_ratio_in, 
#'                                       DEG_list = bulk_DE_cors, case_grep = case_grep,
#'                                       control_grep = control_grep, rda_path = "", 
#'                                       max_proportion_change = 10, print_plots = TRUE, 
#'                                       plot_names = "tst1", theSpecies = "human", 
#'                                       output_directory = "tester",
#'                                       sig_matrix_size = 3000,
#'                                       up_and_downregulated = FALSE, 
#'                                       internet = FALSE)
#' 
#' 
#' @export
#' 
scMappR_and_pathway_analysis <- function(  count_file,signature_matrix, DEG_list, case_grep, control_grep, rda_path = "", max_proportion_change = -9, print_plots=T, plot_names="scMappR",theSpecies = "human", output_directory = "scMappR_analysis",sig_matrix_size = 3000, drop_unknown_celltype = TRUE, internet = TRUE, up_and_downregulated = FALSE, gene_label_size = 0.4, number_genes = -9, toSave=FALSE, newGprofiler = FALSE, path = NULL, deconMethod = "DeconRNASeq") {
  
  
  
  # This function generates cellWeighted_Foldchanges for every cell-type (see deconvolute_and_contextualize) as well as the relative cell-type proportions (which will be reutrned and pushed through)
  # Then, it generates heatmaps of all cellWeighted_Foldchanges, cellWeighted_Foldchanges overlapping with the signature matrix, the signature matrix, the signature matrix overlapping with cellWeighted_Foldchanges
  # Then, if you have WIFI, it will complete g:ProfilR of the reordered cellWeighted_Foldchanges as well as a the ordered list of genes.
  # This function is a wrapper for deconvolute_and_contextualize, as well as pathway_enrich_internal  
  
  # Args:
  # count_file: Normalized RNA-seq count matrix where rows are gene symbols and columns are individuals
  # either the object tself of the path of a TSV file
  # signature_matrix: signature matrix (reccommended odds ratios) of cell-type specificity of genes
  # either the object itself or a pathway to an RData file containing an object named "wilcoxon_rank_mat_or" -- generally internal
  # DEG_list
  # an object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted P-value, and the third the log2FC
  # path to a tsv file containing this info is also acceptable
  #case_grep
  # tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases
  #control_grep
  # tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases
  #max_proportion_change
  # maximum cell-type proportion change -- may be useful if there are many rare cell-types
  # print_plots 
  # whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed
  # The same name of the plots will be completed for top pathways
  #plot_names
  # if plots are being printed, the prefix of their pdf files
  # output_directory
  # the name of the directory that will contain off of the analysis
  #theSpecies
  # -9 if using a precomputed count matrix from scMappR, human otherwise.
  # removes ensembl symbols
  #sig_matrix_size 
  # number of genes in signature matrix for cell-type deconvolution
  # drop_unknown_celltype
  # whether or not to remove "unknown" cell-types from the signature matrix
  # WIFI: Whether you have stable WIFI -- T/F
  # up_and_downregulated: whether you are splitting up/downregulated genes -- T/F
  # gene_label_size = the size of the gene label on the plot
  
  # Returns: 
  # a directory with:
  # cellWeighted_Foldchanges in RData file
  # Cell Type proportions (RData file)
  # cell-type proportions leave one out (RData file)
  # heatmap of cellWeighted_Foldchanges (all)
  # heatmap of cellWeighted_Foldchanges (within signature)
  # heatmap of signature (all)
  # heatmap of signature (overlapping with DEG_list)
  # Pathway enrichment for DEG list(all)
  # RData file and Biological Processes
  # Pathway enrichment of cellWeighted_Foldchanges for each cell-type
  # RData file and biological processes
  
  # load in the count matrix
  
  count_class <- class(count_file)[1] %in% c("character", "data.frame", "matrix")
  
  if(count_class[1] == FALSE) {
    stop("count_file must be of class character, data.frame, or matrix.")
  }
  signature_class <- class(signature_matrix)[1] %in% c("character", "data.frame", "matrix")
  if(signature_class[1] == FALSE) {
    stop("count_file must be of class character, data.frame, or matrix.")
  }
  DEG_list_class <- class(DEG_list)[1] %in% c("character", "data.frame", "matrix") 
  
  if(DEG_list_class[1] == FALSE) {
    stop("DEG_list must be of class character, data.frame, or matrix.")
  }
  case_grep_class <- class(case_grep)[1] %in% c("character", "numeric", "integer")
  case_grep_class[1] == FALSE
  if(case_grep_class[1] == FALSE) {
    stop("case_grep must be of class character (as a single character designating cases in column names) or of class numeric (integer matrix giving indeces of cases).")
  }
  control_grep_class <- class(control_grep)[1] %in% c("character", "numeric", "integer")
  
  if(control_grep_class[1] == FALSE) {
    stop("control_grep must be of class character (as a single character designating controls in column names) or of class numeric (integer matrix giving indeces of controls).")
  }
  
  if(!is.character(theSpecies)) stop("species is not human, mouse, -9, or a species name compatible with g:ProfileR.")
  
  
  if(!is.character(rda_path)) {
    stop("rda_path must be of class list.")
  }
  
  if(!is.numeric(max_proportion_change)) {
    stop("max_proportion_change must be of class numeric.")
  }
  
  if(!is.character(plot_names) ) {
    stop("plot_names must be of class character.")
  }
  
  if(!is.character(output_directory)) {
    stop("output_directory must be of class character.")
  }
  if(!is.numeric(sig_matrix_size) ) {
    stop("sig_matrix_size is not numeric.")
  }
  
  if(!is.numeric(gene_label_size)) {
    stop("gene_label_size must be of class numeric.")
  }
  
  if(!is.numeric(number_genes)) {
    stop("number_genes must be of class numeric.")
  }
  if(all(is.logical(print_plots),is.logical(drop_unknown_celltype),is.logical(internet),is.logical(up_and_downregulated),is.logical(toSave),is.logical(newGprofiler) ) == FALSE) {
    stop("print_plots, drop_unknown_celltype, internet, up_and_down_regulatedtoSave, newGprofiler: must all be class logical.") 
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  theSpecies <- tolower(theSpecies)
  
  if(is.character(count_file)) {
    norm_counts_i <- utils::read.table(count_file, header = TRUE, as.is = TRUE, sep = "\t")
    warning("reading in a count file where the first column is expected to be the row names.")    
    rownames(norm_counts_i) <- norm_counts_i[,1]
    norm_counts_i <- norm_counts_i[,-1]
  } else {
    norm_counts_i <- count_file
  }
  # get the background for later pathway analysis
  background_genes <-  rownames(norm_counts_i)  
  # list of differential expression
  if(is.character(DEG_list)) {
    DEGs <- utils::read.table(DEG_list, header = FALSE, as.is = TRUE, sep = "\t")
  } else {
    DEGs <- as.data.frame(DEG_list)
  }
  
  colnames(DEGs) <- c("gene_name", "padj", "log2fc")
  DEGs$gene_name <- tochr(DEGs$gene_name)
  DEGs$padj <- toNum(DEGs$padj)
  DEGs$log2fc <- toNum(DEGs$log2fc)
  rownames(DEGs) <- DEGs$gene_name
  if(number_genes == -9) {
    number_genes <- as.numeric(nrow(DEGs))
  }
  if((!(is.character(case_grep)) | length(case_grep) > 1)[1]) {
    message("Assuming that case_grep are indeces of 'control'.")
    message("Appending 'scMappR_control to controls.")
    colnames(count_file)[case_grep] <- paste0("scMappR_case_", colnames(count_file)[case_grep])
    case_grep <- "scMappR_case"
  }
  if((!(is.character(control_grep)) | length(control_grep) > 1)[1]) {
    message("Assuming that control_grep indeces of  'control'.")
    message("Appending 'scMappR_control to controls.")
    colnames(count_file)[control_grep] <- paste0("scMappR_control_", colnames(count_file)[control_grep])
    control_grep <- "scMappR_control"
  }
  
  cases <- grep(case_grep, colnames(count_file))
  control <- grep(control_grep, colnames(count_file))
  inter_case_control <- intersect(cases, control)
  if(length(inter_case_control) != 0) {
    stop("Samples in 'case' and 'control' overlap. Please check column names.")
  }
  if(any(length(cases) < 2, length(control) < 2)[1]) {
    stop("There is fewer than two cases or controls, please check 'case_grep' or 'control_grep'.")
  }
  
  if(is.character(signature_matrix)) {
    # assuming that the signature matrix is an RData file containg an object named wilcoxon_rank_mat_or (internal, generally)
    wilcoxon_rank_mat_or <- ""
    load(signature_matrix)
    signature_matrix <- wilcoxon_rank_mat_or
  }
  if(length(grep("-", rownames(signature_matrix))) / length(rownames(signature_matrix)) > 0.75 ) {  
    warning("More than 50 genes contian '-', and the signature matrix is considered internal")
    load(paste0(rda_path,"/bioMart_ortholog_human_mouse.rda"))
    # data(bioMart_ortholog_human_mouse)    
    RN_2 <- get_gene_symbol(signature_matrix)
    
    rownames(signature_matrix) <- RN_2$rowname
    #stop("testing")
    internal_species <- RN_2$species
    if(internal_species != theSpecies) {
      warning(paste0("the species from the signature matrix, ", internal_species, ", does not equal the initially inputted species, ", theSpecies, ". Converting gene symbols of 1-1 orthologs"))
      RN <- rownames(signature_matrix) # gene symbols
      
      thefiles <- list.files(path = rda_path, "bioMart_ortholog_human_mouse.rda")
      
      
      if(length(thefiles) == 0) {
        warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
        #
        datafile <- "bioMart_ortholog_human_mouse.rda"
        metafile <- paste0(datafile)
        url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                      metafile, "?raw=true")
        destfile <- file.path(tempdir(), metafile)
        downloader::download(url, destfile = destfile, mode = "wb")
        load(destfile)
        #
      } else {
        load(paste0(rda_path,"/bioMart_ortholog_human_mouse.rda"))
      }      
      
      # data(bioMart_ortholog_human_mouse) 
      
      rownames(bioMart_orthologs) <- bioMart_orthologs[,internal_species]
      intersected_genes <- intersect(RN, rownames(bioMart_orthologs)) # gene sybmols in signature
      signature_matrix1 <- signature_matrix[intersected_genes,]
      bioMart_orthologs1 <- bioMart_orthologs[intersected_genes,]
      rownames(signature_matrix1) <- bioMart_orthologs1[,theSpecies] # replacing rownames with the species you want
      message("Number of genes: ")
      message(nrow(signature_matrix1))
      message("Number of cell-types: ")
      message(ncol(signature_matrix1))
      signature_matrix <- signature_matrix1
    }
  }
  
  # Identify cell weighted Fold-changes 
  if(toSave == FALSE) {
    print_plots <- FALSE
  }  
  message("Keeping genes in signature matrix that overlap with count matrix.")
  signature_matrix <- signature_matrix[rownames(signature_matrix) %in% rownames(count_file),]
  signature_vars <- apply(signature_matrix, 2, stats::var)
  signature_rowvars <- apply(signature_matrix, 1, stats::var)
  signature_matrix <- signature_matrix[signature_rowvars > 0,]
  signature_matrix <- signature_matrix[,signature_vars > 0]
  message("cell-types with markers that overlap with inputted count matrix")
  message(colnames(signature_matrix))
  
  if((length(unique(colnames(signature_matrix))) < length(colnames(signature_matrix)))[1]) {
    warning("cell-type naming is not unique, appending unique identifier (1:ncol(signature))")
    colnames(signature_matrix) <- paste0(colnames(signature_matrix), "_", 1:ncol(signature_matrix))
  }
  
  cellWeighted_Foldchanges <- deconvolute_and_contextualize(count_file, signature_matrix, DEG_list, case_grep , control_grep, max_proportion_change = max_proportion_change, print_plots = print_plots, plot_names = plot_names, theSpecies = theSpecies, sig_matrix_size = sig_matrix_size, drop_unknown_celltype = drop_unknown_celltype, toSave = toSave, path = path, deconMethod = deconMethod)
  
  # Computing t-test for changes in cell-type proportion
  ttest_decon <- function(x) {
    # This function takes the cell-type proporitons of cases and controls and completes t-tests to see if there are significant changes.  
    Cases <- cellWeighted_Foldchanges$cellType_Proportions[grep(case_grep, rownames(cellWeighted_Foldchanges$cellType_Proportions)),x] # proportion of cell-typeof cases
    Control <- cellWeighted_Foldchanges$cellType_Proportions[grep(control_grep, rownames(cellWeighted_Foldchanges$cellType_Proportions)),x] #controls
    case_control <- stats::t.test(Cases,Control)# t test
    case_control$p.value # p val
    case_control$statistic # t stat
    case_mean <- mean(Cases) # mean cases
    control_mean <- mean(Control) # mean control
    theName <- colnames(cellWeighted_Foldchanges$cellType_Proportions)[x] # cell-type
    toReturn <- c(unname(case_control$p.value), unname(case_control$statistic),case_mean, control_mean )
    return(toReturn)
  } 
  stacked <- lapply(1:ncol(cellWeighted_Foldchanges$cellType_Proportions), ttest_decon)
  proportion_T <- do.call("rbind",stacked) # t test of proportions of each cell-type and staack
  rownames(proportion_T) <- colnames(cellWeighted_Foldchanges$cellType_Proportions)
  colnames(proportion_T) <- c("P.Value", "T.Statistic", "CaseMean", "ControlMean")
  
  cellWeighted_Foldchanges$ProportionT.test <- proportion_T
  
  message(paste0("Making scMappR output directory named", output_directory, "."))
  if(toSave == FALSE) {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended. Returning cellWeighted_Foldchanges and no pathway analysis.")
    message("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended. Returning cellWeighted_Foldchanges and no pathway analysis.")
    return(cellWeighted_Foldchanges)    
  }
  dir.create(paste0(path,"/",output_directory))
  
  scMappR_vals <- cellWeighted_Foldchanges$cellWeighted_Foldchange # scMappR values
  T_test_outs <- cellWeighted_Foldchanges$ProportionT.test
  message("Writing summary of cell-type proportion changes between case and control.")
  utils::write.table(T_test_outs, file = paste0(path,"/",output_directory, "/", plot_names, "_cell_proportion_changes_summary.tsv"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  
  #message(scMappR_vals)
  save(scMappR_vals, file = paste0(path,"/",output_directory, "/",plot_names, "_cellWeighted_Foldchanges.RData"))
  
  cell_proportions_all <- cellWeighted_Foldchanges$cellType_Proportions # all gene CT proportion
  save(cell_proportions_all, file = paste0(path,"/",output_directory, "/",plot_names, "_celltype_proportions.RData"))
  
  leave_one_out_proportions <- cellWeighted_Foldchanges$leave_one_out_proportions # leave one out avg CT proportions
  save(leave_one_out_proportions, file = paste0(path,"/",output_directory, "/",plot_names, "_leaveOneOut_gene_proportions.RData"))
  
  signature_mat <- cellWeighted_Foldchanges$processed_signature_matrix # processed_signaure_matrix
  sigmat_row <- apply(signature_mat, 1, stats::var)
  sigmat_col <- apply(signature_mat, 2, stats::var)
  signature_mat <- signature_mat[which(sigmat_row > 0), which(sigmat_col > 0)]
  save(signature_mat, file = paste0(path,"/",output_directory, "/",plot_names, "_leaveOneOut_gene_proportions.RData"))
  if(nrow(DEG_list) == 1) {
    warning("You only have 1 DEG, no heatmaps can be made. Returning cellWeighted_Foldchange")
    message("You only have 1 DEG, no heatmaps can be made. Returning cellWeighted_Foldchange")
    return(scMappR_vals)
  }
  myheatcol <- grDevices::colorRampPalette(c("lightblue", "white", "orange"))(256)
  
  # generate heatmaps for DEGs
  cex = gene_label_size
  
  scMappR_vals_vars <- apply(scMappR_vals,1,stats::var)
  scMappR_vals <- scMappR_vals[scMappR_vals_vars > 0,]
  
  inter <- intersect(rownames(scMappR_vals),rownames(signature_mat))
  
  #Absolute value plotting of scMappR vals
  if(length(inter) > 3) {
    scMappR_intersect <- scMappR_vals[inter,]
    scMappR_pref <- scMappR_vals[inter,]
    
    
    grDevices::pdf(paste0(path,"/",output_directory,"/",plot_names,"_abs_cwFC_signaure_matrix.pdf")) 
    pl_abs <- pheatmap::pheatmap(abs(as.matrix(scMappR_intersect)), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
    dev.off()
    
    grDevices::pdf(paste0(path,"/",output_directory,"/",plot_names,"abs_cwFC.pdf"))
    pheatmap::pheatmap(abs(as.matrix(scMappR_vals[,pl_abs$tree_col$order])), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10, cluster_cols = FALSE)
    dev.off()
    
    grDevices::pdf(paste0(path,"/",output_directory,"/",plot_names,"_DEG_signature_matrix.pdf")) 
    pheatmap::pheatmap(as.matrix(scMappR_pref[pl_abs$tree_row$order,pl_abs$tree_col$order]), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10, cluster_rows = FALSE, cluster_cols = FALSE)
    dev.off()
    
    
  }
  
  scMappR_vals_up <- as.matrix(scMappR_vals[apply(scMappR_vals,1, sum) > 0,])
  scMappR_vals_down <- as.matrix(scMappR_vals[apply(scMappR_vals,1, sum) < 0,])
  #grDevices::pdf(paste0(path,"/",output_directory,"/",plot_names,"_cell_proportions_heatmap.pdf")) 
  #gplots::heatmap.2(as.matrix(cellWeighted_Foldchanges$cellType_Proportions), Rowv = TRUE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
  #pheatmap::pheatmap(as.matrix(cellWeighted_Foldchanges$cellType_Proportions), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
  #grDevices::dev.off()
  if((nrow(scMappR_vals_up) > 2 & ncol(scMappR_vals_up) > 2)[1]) {
    grDevices::pdf(paste0(path,"/",output_directory, "/", plot_names,"_cwFC_upregulated_DEGs.pdf"))
    
    pheatmap::pheatmap(as.matrix(abs(scMappR_vals_up)), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
    #gplots::heatmap.2(as.matrix(abs(scMappR_vals_up)), Rowv = TRUE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
    grDevices::dev.off()
  } else {
    warning("There were fewer than two upregulated DEGs, therefore a heatmap could not be made.")
    message("There were fewer than two upregulated DEGs, therefore a heatmap could not be made.")
    
  }
  message("Number of genes: ")
  message(nrow(scMappR_vals_down))
  message("Number of cell-types: ")
  message(ncol(scMappR_vals_down))
  if((nrow(scMappR_vals_down) > 2 & ncol(scMappR_vals_down) > 2)[1]) {
    grDevices::pdf(paste0(path,"/",output_directory, "/", plot_names,"_cwFC_downregulated_DEGs.pdf"))
    #gplots::heatmap.2(as.matrix(abs(scMappR_vals_down)), Rowv = TRUE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
    pheatmap::pheatmap(as.matrix(abs(scMappR_vals_down)), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
    
    grDevices::dev.off()
  } else {
    warning("There were fewer than two downregulated DEGs, therefore a heatmap could not be made.")
    message("There were fewer than two downregulated DEGs, therefore a heatmap could not be made.")
    
  }
  
  grDevices::pdf(paste0(path,"/",output_directory, "/",plot_names,"_signature_matrix.pdf"))
  #gplots::heatmap.2(as.matrix(signature_mat), Rowv = TRUE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
  pheatmap::pheatmap(as.matrix(signature_mat), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
  
  grDevices::dev.off()
  
  celltype_preferred_degs <- intersect(rownames(scMappR_vals), rownames(signature_matrix)) # intersect DEGs and Genes in signature matrix
  
  if(length(celltype_preferred_degs) < 3) {
    warning("Fewer than 3 genes are both De and in the signature matrix. Therefore, these heatmaps will not be generated. Furthermore, there is insufficient re-ranking of genes for different pathway analyses to be neccesary. Therefore, here, cellWeighted_Foldchanges are more representative of a scaling factor for each cell-type.")
    message("Fewer than 3 genes are both De and in the signature matrix. Therefore, these heatmaps will not be generated. Furthermore, there is insufficient re-ranking of genes for different pathway analyses to be neccesary. Therefore, here, cellWeighted_Foldchanges are more representative of a scaling factor for each cell-type.")
    
    return(list(cellWeighted_Foldchanges = cellWeighted_Foldchanges))
    
  } else {
    
    # generating the heatmaps for cellWeighted_Foldchanges and signature matrix odds ratios that overlap with one another
    up_signature <- intersect(rownames(scMappR_vals_up), celltype_preferred_degs)
    down_signature <- intersect(rownames(scMappR_vals_down), celltype_preferred_degs)
    
    
    
    
    
    signature_mat_up <- as.matrix(signature_mat[up_signature,])
    signature_mat_down <- as.matrix(signature_mat[down_signature,])
    scMappR_vals_up1 <- as.matrix(scMappR_vals_up[up_signature,])
    scMappR_vals_down1 <- as.matrix(scMappR_vals_down[down_signature,])
    
    # Upregulated DEG Heatmap
    if(nrow(signature_mat_up) > 2 & ncol(signature_mat_up) > 2) {
      grDevices::pdf(paste0(path,"/",output_directory, "/",plot_names,"cwFC_signature_matrix_upregulated_DEGs.pdf"))
      #pl <- gplots::heatmap.2(as.matrix(signature_mat_up), Rowv = TRUE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
      pl <- pheatmap::pheatmap(as.matrix(scMappR_vals_up1), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
      
      #message(pl)
      grDevices::dev.off()
      
      
      grDevices::pdf(paste0(path,"/",output_directory, "/", plot_names,"_signature_matrix_upregulated_DEGs.pdf"))
      #gplots::heatmap.2(as.matrix(scMappR_vals_up1[rev(colnames(pl$carpet)),pl$colInd]),Colv=F, Rowv = FALSE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
      tst <- pheatmap::pheatmap(as.matrix(signature_mat_up)[pl$tree_row$order,pl$tree_col$order], color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10, cluster_cols = FALSE, cluster_rows = FALSE)
      
      grDevices::dev.off()
    } else {
      warning("There were fewer than two cell-type specific, upregulated DEGs, therefore a heatmap could not be made.")
      message("There were fewer than two cell-type specific, upregulated DEGs, therefore a heatmap could not be made.")
      
    }
    
    # Downregulated DEG Heatmap
    
    if((nrow(signature_mat_down) > 2 & ncol(signature_mat_down) > 2)[1]) {
      grDevices::pdf(paste0(path,"/",output_directory, "/",plot_names,"_cwFC_signature_matrix_downregulated_DEGs.pdf"))
      pl2 <- pheatmap::pheatmap(as.matrix(scMappR_vals_down1), color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10)
      grDevices::dev.off()
      
      
      grDevices::pdf(paste0(path,"/",output_directory, "/", plot_names,"_signature_matrix_downregulated_DEGs.pdf"))
      #gplots::heatmap.2(as.matrix(abs(scMappR_vals_down1[rev(colnames(pl2$carpet)),pl2$colInd])),Colv=F, Rowv = FALSE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
      tst <- pheatmap::pheatmap(as.matrix(signature_mat_down)[pl2$tree_row$order,pl2$tree_col$order], color = myheatcol, scale = "row", fontsize_row = cex, fontsize_col = 10, cluster_cols = FALSE, cluster_rows = FALSE)
      
      grDevices::dev.off()
    }  else {
      warning("There were fewer than two cell-type specific, downregulated DEGs, therefore a heatmap could not be made.")
      message("There were fewer than two cell-type specific, downregulated DEGs, therefore a heatmap could not be made.")
      
    }
  }
  if(internet == FALSE) {
    warning("There is not a reported stable internet (WIFI = FALSE) and therefore pathway analysis with g:Prof")
    return("Done!")
  }
  
  up_and_down_together <- pathway_enrich_internal(  DEGs, theSpecies, scMappR_vals, background_genes, output_directory, plot_names, number_genes = number_genes, toSave=TRUE,path=path, newGprofiler = newGprofiler)
  if(up_and_downregulated == TRUE)  {
    message("Splitting genes by up- and down-regulated and then repeating analysis")
    rownames(DEGs) <- DEGs$gene_name
    upGenes <- DEGs$gene_name[DEGs$log2fc > 0]
    downGenes <- DEGs$gene_name[DEGs$log2fc < 0]
    
    
    
    upDEGs <- DEGs[upGenes,]
    downDEGs <- DEGs[downGenes,]
    
    if(length(upGenes) > 2) {
      upDir <- paste0(output_directory,"/upregulated") 
      dir.create(paste0(path,"/",upDir))
      message("Pathway analysis of upregulated genes")
      upcellWeighted_Foldchanges <- scMappR_vals[upGenes,]
      up_only <- pathway_enrich_internal(  upDEGs, theSpecies, upcellWeighted_Foldchanges, background_genes, upDir, plot_names, number_genes = number_genes, toSave=TRUE, path = path, newGprofiler = newGprofiler)
    } else {
      warning("There are fewer than three upregulated DEGs.")
    }
    message("Pathway analysis of downregulated genes")
    if(length(downGenes) > 2) {
      downDir <- paste0(output_directory, "/downregulated")
      
      dir.create(paste0(path,"/",downDir))
      DowncellWeighted_Foldchanges <- scMappR_vals[downGenes,]
      down_only <- pathway_enrich_internal(  downDEGs, theSpecies, DowncellWeighted_Foldchanges, background_genes, downDir, plot_names, number_genes = number_genes, toSave=TRUE, path = path, newGprofiler = newGprofiler)    
    } else {
      warning("There are fewer than three downregulated DEGs.")
      
    }
  }
  
  return(list(cellWeighted_Foldchanges = cellWeighted_Foldchanges, paths = up_and_down_together))
}
