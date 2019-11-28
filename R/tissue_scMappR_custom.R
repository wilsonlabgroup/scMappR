#' Gene List Visualization and Enrichment with Custom Signature Matrix
#'
#' This function visualizes signature matrix, cluster's subsetted genes, completes enrichment of individual cell-types and co-enrichment.
#' 
#' This function is roughly the same as tissue_scMappR_internal, however now there is a custom signature matrix.
#' it generates a heatmap of the signature matrixand your inputted gene list as well as single cell-type and 
#' co-celltype enrichment.
#'
#'
#' @rdname tissue_scMappR_custom
#' @name tissue_scMappR_custom
#'
#' @param gene_list A list of gene symbols, mouse or human.
#' @param signature_matrix precomputed signature matrix with matching gene names.
#' @param output_directory Directory made containing outputs
#' @param gene_cutoff value cutoff (generally rank := log10(Padj)) for a gene to be considered a marker
#' @param is_pvalue If signature matrix is p-value before rank is applied (not recommended ) (T/F)
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#'
#' @return \code{tissue_scMappR_custom} A list containing the entire signature matrix, the matrix subsetted for your genes, enrichment of each cell-type, and co-enrichment. \cr
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
#'
#' @examples 
#' \notrun {
#' 
#' # load in signature matrices
#' load("data/Preoptic_region_example.rda")
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' # data(Preoptic_region_example)
#' Signature <- POA_Rank_signature
#'  rowname <- get_gene_symbol(Signature)
#'  rownames(Signature) <- rowname$rowname
#'  genes <- rownames(Signature)[1:200]
#'  
#'  # Assuming Signature_matrices_pVal.rda is in you "~/Documents/scMappR/data" directory
#'  # rda_path1 = "~/Documents/scMappR/data"
#'  internal <- tissue_scMappR_custom(genes,Signature,output_directory = "scMappR_Test_custom")
#'  }
NULL
#' @rdname tissue_scMappR_custom
#' @export
#' 
tissue_scMappR_custom <- function(gene_list, signature_matrix ,output_directory = "custom_test", toSave = FALSE, gene_cutoff = 1, is_pvalue = TRUE) {
  # This function is roughly the same as tissue_scMappR_internal, however now there is a custom signature matrix.
  # it generates a heatmap of the signature matrixand your inputted gene list as well as single cell-type and 
  # co-celltype enrichment.
  
  # Args:
  # a list of symbols that matches their background
  # a background signature matrix
  # output_directory: the name of the directory it will generate
  # gene cutoff: the cutoff for a gene to be includedin the cell-type
  # is_pvalue: if the value cutoff is a p-value, we will apply a -log10() transformatiion
  
  # Returns: generates a directory containing
  # heatmap for signature matrix and matrix intersecting with heatmap
  # RData file containing cell-type preferences
  # tsv files with gene set enrichment for single cell-types and co-enrichment
  # It also returns the object of the single cell-type preferences
  
  
  outDir <- paste0(output_directory)
  study_names <- outDir
  if(toSave == TRUE) {
    dir.create(outDir)
  } else {
    warning("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.")
    print("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.", quote = F)
    
  }
  single_cell_studies <- list()
  i=1 # makes the gene set enrichment syntax compatible with heatmap_generation_internal
  background = signature_matrix
  background_genes <- rownames(signature_matrix)
  
  background_heatmap <- heatmap_generation(background_genes, comp = paste0(outDir, "/", study_names,"_background"), reference = background, isBackground = TRUE, pVal = gene_cutoff, isPval = is_pvalue, toSave = toSave)  
  # heatmap generation of the entire signature matrix
  gene_list_heatmap <- heatmap_generation(gene_list, comp = paste0(outDir, "/", study_names,"_genelist"), reference = background, pVal = gene_cutoff, isPval = is_pvalue, toSave = toSave)
  if(class(gene_list_heatmap) == "character") {
    warning("0 or 1 input genes were cell-type specific. No downstream analysis available.")
    warning("With 0 genes being cell-type specific, I would make sure that they are the same gene symbols.")
    return("0 or 1 input genes were cell-type specific. No downstream analysis available.")
  }
  # heatmap generation of the background matrix as well as getting preferred genes based on your cutoff, p-value or otherwise
  singleCTpreferences <- single_gene_preferences(gene_list_heatmap, background_heatmap, study_names, outDir = output_directory, toSave = toSave)
  
  if(toSave == TRUE) {
    write.table(singleCTpreferences, file = paste0(outDir,"/",outDir, "_celltype_preferences.tsv"), quote=F, row.names = F, col.names = T, sep = "\t")
  }
  # cell-type preferences for indidual cell-types
  sig <- singleCTpreferences[singleCTpreferences$pFDR < 0.05,]
  sig$Odds_Ratio <- toNum(sig$Odds_Ratio)
  sig <- sig[sig$Odds_Ratio > 1,]
  
  if(nrow(sig) <2 ) {
    # if fewer than two cell-types are enriched
    print("co-enrichment cannot be measured as one or fewer CTs are enriched")
    coCTpreferences <- "co-enrichment cannot be measured as one or fewer CTs are enriched"
    output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
    single_cell_studies <- output 
    
  } else {
    sig <- sig[order(sig$pFDR),]
    
    coCTpreferences <- coEnrich(sig, gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory, toSave = toSave)
    # test to see if the same genes are responsible for the enrichment of multiple cell-types
    output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
    
    single_cell_studies <- output 
  }
  names(single_cell_studies) <- c("background_heatmap", "gene_list_heatmap","single_celltype_preferences", "group_celltype_preferences")
  return(single_cell_studies)
}

