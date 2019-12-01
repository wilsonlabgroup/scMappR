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
#' @param species Human, mouse, or a charcter that is compatible with gProfileR.
#' @param background A list of background genes to test against.
#' @param output_directory Path to the directory where files will be saved.
#' @param plot_names Names of output.
#' @param number_genes Number of genes to if there are many many DEGs
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#' 
#' @return \code{pathway_enrich_internal} Plots and pathway enrichment of bulk DE and STVs. \cr
#' 
#' @import matrixStats
#' @import DeconRNASeq
#' @import S4Vectors
#' @import ggplot2
#' @import gplots
#' @import graphics
#' @import Seurat
#' @import GSVA
#' @import stats
#' @import utils
#' @import downloader
#' @import grDevices
#'
#' @examples 
#' \donttest {
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
#' norm <- deconvolute_and_contextualize(bulk_normalized, odds_ratio_in, bulk_DE_cors, case_grep = case_grep,
#'                                       control_grep = control_grep, max_proportion_change = max_proportion_change,
#'                                       print_plots = print_plots, theSpecies = theSpecies)
#' background = rownames(bulk_normalized)
#' dir.create("test_path")
#' pathway_enrich_internal(bulk_DE_cors, "human", norm$scMappR_transformed_values,
#'                         background, "test_path", "test_figs", toSave = TRUE)
#' 
#' }
#' @export
#' 
pathway_enrich_internal <- function(DEGs, theSpecies, scMappR_vals, background_genes, output_directory, plot_names, number_genes, toSave = FALSE) {
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
  
  print("Reordering DEGs from bulk dataset.", quote = F)
  DEG_Names <- rownames(DEGs)[order(DEGs$padj)]
  if(theSpecies == "human") species_bulk <- "hsapiens"
  if(theSpecies == "mouse") species_bulk <- "mmusculus"
  
  ### We could theoretically split this into it's own function if we wanted
  # Pathway enrichment 
  ordered_back_all <- gProfileR::gprofiler(DEG_Names, species_bulk, ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
  ordered_back_all_tf <- gProfileR::gprofiler(DEG_Names, species_bulk, ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")

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
  
  print("Compelting pathway analysis of STVs.", quote = F)
  paths <- gProfiler_STV(scMappR_vals,species =  theSpecies, background = background_genes, gene_cut = number_genes) # re-rodered pathway analysis
  
  biological_pathways <- paths$BP # Biological pathways
  transcription_factors <- paths$TF # Transcription factors
  save(biological_pathways, file = paste0(output_directory, "/",plot_names,"_reordered_pathways.RData"))
  save(transcription_factors, file = paste0(output_directory, "/",plot_names,"_reordered_transcription_factors.RData"))
  
  print("Plotting top 10 biological proccesses and TFs", quote= F)
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
