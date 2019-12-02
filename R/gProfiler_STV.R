#' Pathway enrichment for STV's
#' 
#' This function runs through each list of STV's and completes both pathway and TF enrichment.
#'
#' This function takes a matrix of scMappR_Transformed_Values and a species (Human, mouse, or a character directly compatible with g:ProfileR).
#' Before completing pathway analysis with g:ProfileR. Enriched pathways are stored in a list and returned.
#'
#' @rdname gProfiler_STV
#' @name gProfiler_STV
#'
#' @param STV_matrix Matrix of scMappR Transformed Values from the deconvolute_and_contextualize functions.
#' @param species Human, mouse, or a charcter that is compatible with gProfileR.
#' @param background A list of background genes to test against.
#' @param gene_cut The top number of genes in pathway analysis.
#' 
#' @return \code{gProfiler_STV} A List of significantly enriched pathways and TFs (correction_method = FDR, hier_sorting = moderate), for every cell-type. \cr
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
#' @import gProfileR
#'
#' @examples 
#' \donttest{
#' 
#' data(PBMC_scMappR)
#' 
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' 
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' norm <- deconvolute_and_contextualize(bulk_normalized, odds_ratio_in,
#'                                       bulk_DE_cors, case_grep = case_grep, 
#'                                       control_grep = control_grep,
#'                                        max_proportion_change = max_proportion_change,
#'                                        print_plots = print_plots,
#'                                        theSpecies = theSpecies)
#' background = rownames(bulk_normalized)
#' STVs <- gProfiler_STV(norm$scMappR_transformed_values, theSpecies, background)
#' 
#'  }
#'  
#' @export
#' 



gProfiler_STV <- function(STV_matrix, species , background , gene_cut ) {
  
  # Args: 
  # STV_matrix: matrix of scMappR Transformed Values from the deconvolute_and_contextualize functions
  # species: human, mouse, or a charcter that is compatible with gProfileR
  # background: a list of background genes to test against
  # Returns:
  # A List of significantly enriched pathways and TFs (correction_method = FDR, hier_sorting = moderate), for every cell-type
  
  
  if(species == "human") {
    print("Assuming species = Human", quote = F)
    theSpecies <- "hsapiens"
    
  }
  if(species == "mouse") {
    print("Assuming species = Mouse", quote = F)
    theSpecies <- "mmusculus"
    
    
  }
  if(!(species %in% c("human", "mouse"))) {
    warning("We are assuming that you inputted 'species' as something directly compatible with g:ProfileR") 
    theSpecies <- species
    
  }
  gProfiler_internal <- function(x, theSpecies1 = theSpecies, background_genes = background, cutgenes = gene_cut) {
    #Internal -- This function takes a signle STV and computes gProfiler BP and TF enrichment with the list re-ordered for each cell-type
    # Args:
    # X: a numeric index giving the cell-type in the STV matrix
    # theSpecies: the species that will be used for gProfiler
    # background_genes: the detected genes in your RNA-seq experiment (what you want to use as a background)
    # Returns:
    # significantly enriched BPs or TFs for at least one cell-type
    
    
    print(paste0("Re-ordering by absolute value of STVs on cell-type ", colnames(STV_matrix)[x]))
    STV_matrix1 <- STV_matrix[order(abs(STV_matrix[,x]), decreasing = T) ,]
    
    STV_matrix1 <- rownames(STV_matrix1)[abs(STV_matrix1[,x]) > 1e-10]
    if(cutgenes != -9) {
      print(paste0("Taking the top ", cutgenes, " genes for pathway analysis."), quote = F)
      STV_matrix1 <- STV_matrix1[1:cutgenes]
    }
    ordered_back_all <- gProfileR::gprofiler(STV_matrix1, theSpecies1, ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
    ordered_back_all_tf <- gProfileR::gprofiler(STV_matrix1, theSpecies1, ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = background_genes, correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
    return(list(BPs = ordered_back_all, TFs = ordered_back_all_tf))
    
  }
  print(colnames(STV_matrix))
  print(theSpecies)
  paths <- lapply(1:ncol(STV_matrix), gProfiler_internal)
  
  BPs <- TFs <- list()
  # Converting BPs and TFs for every cell-type into a list
  for(i in 1:ncol(STV_matrix)) {
    BPs[[i]] <- paths[[i]]$BPs
    TFs[[i]] <- paths[[i]]$TFs
  }
  names(BPs) <- names(TFs) <- colnames(STV_matrix)
  return(list(BP = BPs, TF = TFs))
  
}
