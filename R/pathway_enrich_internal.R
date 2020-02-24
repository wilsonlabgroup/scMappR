#' Pathway enrichment for cellWeighted_Foldchange's and bulk gene list
#' 
#' This function completes pathway enrichment of cellWeighted_Foldchanges and bulk gene list.
#'
#' Internal: Pathway analysis of differentially expressed genes (DEGs) and cell weighted Fold-changes (cellWeighted_Foldchanges) for each cell-type. Returns .RData objects of differential analysis as well as plots of the top bulk pathways.
#' It is a wrapper for making barplots, bulk pathway analysis, and gProfiler_cellWeighted_Foldchange.
#' 
#' @rdname pathway_enrich_internal
#' @name pathway_enrich_internal
#' 
#' @param DEGs Differentially expressed genes (gene_name, padj, log2fc).
#' @param theSpecies Human, mouse, or a charcter that is compatible with gProfileR.
#' @param scMappR_vals cell weighted Fold-changes of differentially expressed genes.
#' @param background_genes A list of background genes to test against.
#' @param output_directory Path to the directory where files will be saved.
#' @param plot_names Names of output.
#' @param number_genes Number of genes to if there are many, many DEGs.
#' @param newGprofiler Whether to use gProfileR or gprofiler2 (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#' 
#' @return \code{pathway_enrich_internal} Plots and pathway enrichment of bulk DEGs and cellWeighted_Foldchanges. \cr
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#'
#' @examples 
#' \donttest{
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
#' toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
#'                                       bulk_DE_cors, case_grep = case_grep,
#'                                       control_grep = control_grep, rda_path = "", 
#'                                       max_proportion_change = 10, print_plots = TRUE, 
#'                                        plot_names = "tst1", theSpecies = "human", 
#'                                        output_directory = "tester",
#'                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
#'                                        internet = FALSE)
#' 
#' }
#' @export
#' 
pathway_enrich_internal <- function(DEGs, theSpecies, scMappR_vals, background_genes, output_directory, plot_names, number_genes = -9,  newGprofiler = FALSE, toSave = FALSE, path = NULL) {
  # Internal: Pathway analysis od DEGs and cellWeighted_Foldchanges for each cell-type. Returns RData objects of differential analysis as well as plots of the top bulk pathways.
  # It is a wrapper for making barplots, bulk pathway analysis, and gProfiler_cellWeighted_Foldchange
  # Args:
  # DEGs = differentially expressed genes
  # theSpecies: human/mouse or something compatible with gProfileR
  # scMappR_vals: matrix of of cellWeighted_Foldchanges for each cell-type.
  # background genes = the genes in the count dataset
  # output directory = path to the directory where files will be saved
  # plot_names = names of output
  # Returns:
  # plots and pathway enrichment of bulk DE and cellWeighted_Foldchanges.
  
  DEGs_class <- class(DEGs)[1] %in% c("data.frame", "matrix")
  if(DEGs_class[1] == FALSE) {
    stop("DEGs must be of class data frame or matrix.")
  }
  
  scMappR_vals_class <- class(scMappR_vals)[1] %in% c("data.frame", "matrix")
  if(scMappR_vals_class[1] == FALSE) {
    stop("scMappR_vals must be a data.frame or matrix.")
  }
  if(!is.character(background_genes)) {
    stop("background_genes must be a character vector of gene names.")
  }
    
  if(!is.character(theSpecies)) {
    stop("the species must be a character, human, mouse, or a species compatible with gprofiler.")
  }
  
  if(!is.character(output_directory)) {
    stop("output_directory must be a character, human, mouse, or a species compatible with gprofiler.")
  }
  
  if(!is.character(plot_names)) {
    stop("plot_names must be a character, human, mouse, or a species compatible with gprofiler.")
  }
  if(!is.numeric(number_genes)) {
    message(number_genes)
    message(class(number_genes))
    stop("number_genes must be of class numeric.")
  }
  if(!is.logical(newGprofiler)) {
    stop("newGprofiler must be logical TRUE/FALSE.")
  }
  
  
  if(toSave == FALSE) {
    stop("toSave = FALSE and therefore scMappR is not allowed to print pathways. For this function to work, please set toSave = TRUE")
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  message("Reordering DEGs from bulk dataset.")
  DEG_Names <- rownames(DEGs)[order(DEGs$padj)]
  if(theSpecies == "human") species_bulk <- "hsapiens"
  if(theSpecies == "mouse") species_bulk <- "mmusculus"
  
  ### We could theoretically split this into it's own function if we wanted
  # Pathway enrichment 

  if(newGprofiler == TRUE) {
  
  ordered_back_all <- gprofiler2::gost(query = DEG_Names, organism = species_bulk, ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "fdr",  custom_bg =background_genes, numeric_ns = "", sources = c("GO:BP", "KEGG", "REAC"))  

  if(is.null(ordered_back_all)) { # grpofiler2 returns null if nothing is significant, making compatile with plotBP and make_TF_barplot
  ordered_back_all <- as.data.frame(matrix(c(1," "),nrow = 1))
  colnames(ordered_back_all) <- c("p_value", "term_name")
  ordered_back_all <- ordered_back_all[-1,]
  } else {
    ordered_back_all <- ordered_back_all$result
    ordered_back_all <- ordered_back_all[ordered_back_all$term_size > 15 & ordered_back_all$term_size < 2000 & ordered_back_all$intersection_size > 2,]
  }  

  ordered_back_all_tf <- gprofiler2::gost(query = DEG_Names, organism = species_bulk, ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "fdr", custom_bg =background_genes, numeric_ns = "", sources = c("TF"))  
  
  if(is.null(ordered_back_all_tf)) { # grpofiler2 returns null if nothing is significant, making compatile with plotBP and make_TF_barplot
    ordered_back_all_tf <- as.data.frame(matrix(c(1," "),nrow = 1))
    colnames(ordered_back_all_tf) <- c("p_value", "term_name")
    ordered_back_all_tf <- ordered_back_all_tf[-1,]
  } else {
    ordered_back_all_tf <- ordered_back_all_tf$result
    ordered_back_all_tf <- ordered_back_all_tf[ordered_back_all_tf$term_size > 15 & ordered_back_all_tf$term_size < 5000 & ordered_back_all_tf$intersection_size > 2,]
  }  
  } else {
    ordered_back_all <- gProfileR::gprofiler(DEG_Names, species_bulk, ordered_query = TRUE, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
    ordered_back_all_tf <- gProfileR::gprofiler(DEG_Names, species_bulk, ordered_query = TRUE, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
    ordered_back_all$term_name <- ordered_back_all$term.name
    ordered_back_all$p_value <- ordered_back_all$p.value
    ordered_back_all_tf$term_name <- ordered_back_all_tf$term.name
    ordered_back_all_tf$p_value <- ordered_back_all_tf$p.value
  }
  
  save(ordered_back_all, file = paste0(path, "/", output_directory,"/Bulk_pathway_enrichment.RData"))
  save(ordered_back_all_tf, file = paste0(path, "/", output_directory,"/Bulk_TF_enrichment.RData"))
  #plotting paths
  grDevices::png(file = paste0(path, "/", output_directory,"/Bulk_pathway_enrichment.png"))
  bulk_bp <- plotBP(ordered_back_all)
  message(bulk_bp)
  grDevices::dev.off()
  
  #plotting TFs
  grDevices::png(file = paste0(path,"/",output_directory,"/Bulk_TF_enrichment.png"))
  bulk_bp <- make_TF_barplot(ordered_back_all_tf, top_tf = 10)
  message(bulk_bp)
  grDevices::dev.off()

  
  save(ordered_back_all, file = paste0(path,"/",output_directory, "/",plot_names,"_bulk_pathways.RData"))
  save(ordered_back_all_tf, file = paste0(path,"/",output_directory, "/",plot_names,"_bulk_transcription_factors.RData"))
  
  message("Compelting pathway analysis of cellWeighted_Foldchanges.")
  paths <- gProfiler_cellWeighted_Foldchange(scMappR_vals,species =  theSpecies, background = background_genes, gene_cut = number_genes, newGprofiler = newGprofiler) # re-rodered pathway analysis
  
  biological_pathways <- paths$BP # Biological pathways
  transcription_factors <- paths$TF # Transcription factors
  save(biological_pathways, file = paste0(path, "/", output_directory, "/",plot_names,"_reordered_pathways.RData"))
  save(transcription_factors, file = paste0(path, "/", output_directory, "/",plot_names,"_reordered_transcription_factors.RData"))
  
  message("Plotting top 10 biological proccesses and TFs")
  BP_dir <- paste0(path, "/", output_directory, "/BP_barplot")
  TF_dir <- paste0(path, "/", output_directory, "/TF_barplot")
  
  dir.create(BP_dir)
  for(i in 1:length(biological_pathways)) {
    # printing top 10 pathways for each celltype
    grDevices::pdf(file = paste0(BP_dir,"/",plot_names,"_",names(biological_pathways)[i],"_BP.pdf"))
    BP <- plotBP(biological_pathways[[i]], top_bp = 10)
    print(BP)
    grDevices::dev.off()
  }
  dir.create(TF_dir)
  for(i in 1:length(transcription_factors)) {
    # printing top 10 transcriptionfactors for ech celltype
    grDevices::pdf(file = paste0(TF_dir,"/",plot_names,"_",names(transcription_factors)[i],"_TF.pdf"))
    TF  <- make_TF_barplot(transcription_factors[[i]], top_tf = 10)
    print(TF)
    grDevices::dev.off()
  }
  return(list(BPs = biological_pathways, TFs = transcription_factors))
}
