#' Pathway enrichment for STV's and bulk gene list
#' 
#' This function completes pathway enrichment of STVs and bulk gene list.
#'
#' Internal: Pathway analysis od DEGs and STVs for each cell-type. Returns RData objects of differential analysis as well as plots of the top bulk pathways.
#' It is a wrapper for making barplots, bulk pathway analysis, and gProfiler_STV
#' 
#' @rdname pathway_enrich_internal
#' @name pathway_enrich_internal
#' 
#' @param DEGs Differentially expressed genes (gene_name, padj, log2fc).
#' @param theSpecies Human, mouse, or a charcter that is compatible with gProfileR.
#' @param scMappR_vals scMappR Transformed Values of differentially expressed genes.
#' @param background_genes A list of background genes to test against.
#' @param output_directory Path to the directory where files will be saved.
#' @param plot_names Names of output.
#' @param number_genes Number of genes to if there are many many DEGs.
#' @param newGprofiler Whether to use g:ProfileR or gprofiler2 (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' 
#' @return \code{pathway_enrich_internal} Plots and pathway enrichment of bulk DE and STVs. \cr
#' 
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom gplots heatmap.2
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
#' data(PBMC_scMappR)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' norm <- deconvolute_and_contextualize(bulk_normalized, odds_ratio_in, bulk_DE_cors,
#'                                        case_grep = case_grep, control_grep = control_grep,
#'                                         max_proportion_change = max_proportion_change,
#'                                       print_plots = print_plots, theSpecies = theSpecies)
#' background = rownames(bulk_normalized)
#' dir.create("test_path")
#' pathway_enrich_internal(bulk_DE_cors, "human", norm$scMappR_transformed_values,
#'                         background, "test_path", "test_figs", toSave = TRUE)
#' 
#' }
#' @export
#' 
pathway_enrich_internal <- function(DEGs, theSpecies, scMappR_vals, background_genes, output_directory, plot_names, number_genes,  newGprofiler, toSave = FALSE) {
  # Internal: Pathway analysis od DEGs and STVs for each cell-type. Returns RData objects of differential analysis as well as plots of the top bulk pathways.
  # It is a wrapper for making barplots, bulk pathway analysis, and gProfiler_STV
  # Args:
  # DEGs = differentially expressed genes
  # theSpecies: human/mouse or something compatible with g:ProfileR
  # scMappR_vals: matrix of of STVs for each cell-type.
  # background genes = the genes in the count dataset
  # output directory = path to the directory where files will be saved
  # plot_names = names of output
  # Returns:
  # plots and pathway enrichment of bulk DE and STVs.
  
  if(toSave == FALSE) {
    stop("toSave = FALSE and therefore scMappR is not allowed to print pathways. For this function to work, please set toSave = TRUE")
  }
  
  print("Reordering DEGs from bulk dataset.", quote = FALSE)
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
  
  save(ordered_back_all, file = paste0(output_directory,"/Bulk_pathway_enrichment.RData"))
  save(ordered_back_all_tf, file = paste0(output_directory,"/Bulk_TF_enrichment.RData"))
  #plotting paths
  grDevices::pdf(file = paste0(output_directory,"/Bulk_pathway_enrichment.pdf"))
  bulk_bp <- plotBP(ordered_back_all)
  print(bulk_bp)
  grDevices::dev.off()
  
  #plotting TFs
  grDevices::pdf(file = paste0(output_directory,"/Bulk_TF_enrichment.pdf"))
  bulk_bp <- make_TF_barplot(ordered_back_all_tf, top_tf = 10)
  print(bulk_bp)
  grDevices::dev.off()

  
  save(ordered_back_all, file = paste0(output_directory, "/",plot_names,"_bulk_pathways.RData"))
  save(ordered_back_all_tf, file = paste0(output_directory, "/",plot_names,"_bulk_transcription_factors.RData"))
  
  print("Compelting pathway analysis of STVs.", quote = FALSE)
  paths <- gProfiler_STV(scMappR_vals,species =  theSpecies, background = background_genes, gene_cut = number_genes, newGprofiler = newGprofiler) # re-rodered pathway analysis
  
  biological_pathways <- paths$BP # Biological pathways
  transcription_factors <- paths$TF # Transcription factors
  save(biological_pathways, file = paste0(output_directory, "/",plot_names,"_reordered_pathways.RData"))
  save(transcription_factors, file = paste0(output_directory, "/",plot_names,"_reordered_transcription_factors.RData"))
  
  print("Plotting top 10 biological proccesses and TFs", quote= FALSE)
  BP_dir <- paste0(output_directory, "/BP_barplot")
  TF_dir <- paste0(output_directory, "/TF_barplot")
  
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
