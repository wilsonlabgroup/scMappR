#' Generate signature matrix 
#'
#' Convert a list of cell-type markers from FindMarkers in Seurat to a signature matrix defined by odds ratio and rank.
#'
#' Take an list of compiled DEGs from different cell types, identify what the cell-types are using the fisher's exact method, and then convert into a signature matrix for both the adjusted p-value and odds ratio
#' 
#' @rdname generes_to_heatmap
#' @name generes_to_heatmap
#'
#' @param generes A list of cell-tpe markers with fold-changes and P-vlaues (FindMarkers output in Seurat)
#' @param species The species of gene symbols, if not internal, "human" or "mouse"
#' @param naming_preference Likely cell-types given tissues (to be passed into human_mouse_ct_marker_enrich)
#' @param make_names Identify names of cell-type markers using the Fisher's Exact-Test method (T/F).
#' @param internal If this function is pre-processing from panglao (T/F).
#' 
#'  
#'
#' @return \code{generes_to_heatmap} A list containing a signature matrix by rank := -1*log10(Pfdr) and by fold-change (only increasing). Additionally it returns the top (up to) 30 CT markers for each cell-type, as well as the name of each cell-type (from the signature methods method). \cr
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
#'  
#'  
#'  # load generes object
#' load("~/scMappR/data/Preoptic_region_example.rda")
#' #data(Preoptic_region_example)
#'  POA_generes <- POA_example$POA_generes
#'  POA_OR_signature <- POA_example$POA_OR_signature
#'  POA_Rank_signature <- POA_example$POA_Rank_signature
#' signature <- generes_to_heatmap(POA_generes,species = -9, make_names = F)
NULL
#' @rdname generes_to_heatmap
#' @export

generes_to_heatmap <- function(generes,  
                               species = "human",  
                               naming_preference = -9, 
                               rda_path = "~/scMappR/data", 
                               make_names = T, 
                               internal = FALSE) {
  # take an list of compiled DEGs from different cell types, identify what the cell-types are using the fisher's exact method, and then convert into a signature matrix for both the adjusted p-value and odds ratio
  # Args:
  # generes: A list of cell-type markers where each element of the list has the gene symbol, p_adj, and fold-change: the output of the FindMarkers function in Seurat
  # colnames are asfollows
  #  p_val_adj avg_logFC
  # species: human, mouse, or -9 if internal
  # naming_preference: if you have an idea of what cell-types should be in the dataset, then using one of the arguments in the 
  # get_naming_preference() funciton will increase the rank that cell-types falling into that category
  # make_names: If you ant the cell-types to be named or not
  # internal: If this is part of our internal pipeline, we apend the mouse and human symbol onto the final matrix so that we can easily tell afterwards
  # Returns:
  # A list containing a signature matrix by rank := -1*log10(Pfdr) and by fold-change (only increasing). 
  # additionally it returns the top (up to) 30 CT markers for each cell-type, as well as the name of each cell-type (from the signature methods method)
  if(species == -9) {
    # if it's internal and symbols and ENSBML are attached
    for(i in 1:length(generes)) {
      names1 <- get_gene_symbol(generes[[i]])
      rownames(generes[[i]]) <- names1$rowname
    }
    species <- names1$species
  }
  topGenes <- topgenes_extract(generes) # take the top 30 genes
  if(make_names == T) {
    cell_name <- human_mouse_ct_marker_enrich(topGenes, theSpecies = species, cell_marker_path = rda_path) 
    # get the names of each of the cell types
    cell_name <- cell_name$cellTypes # attach the appropriate cell names
  } else {
    cell_name = colnames(generes)
  }
  genes <- c()
  for(i in 1:length(generes)) {
    genes <- c(genes,rownames(generes[[i]]))
  }
  # generate matrix with ranks or odds ratio
  # only genes that are preferentially expressed in at min one cell type
  # is removed
  genes_uni <- unique(genes)
  if(internal == TRUE) {
    if(species == "human") {
      genes_uni <- paste0(genes_uni, "-ENSG00000000000.0")
      for(i in 1:length(generes)) {
        rownames(generes[[i]]) <- paste0(rownames(generes[[i]]),"-ENSG00000000000.0")
      }
      
    }
    if(species == "mouse") {
      genes_uni <- paste0(genes_uni, "-ENSMUSG00000000000.0")
      for(i in 1:length(generes)) {
        rownames(generes[[i]]) <- paste0(rownames(generes[[i]]),"-ENSMUSG00000000000.0")
      }
    }
  }
  
  # generating the signature matrix developed by the rank of p-vlaues
  scmappr <- matrix(0, nrow = length(genes_uni), ncol = length(names(generes))) 
  rownames(scmappr) <- genes_uni
  colnames(scmappr) <- names(generes)
  for(i in 1:length(generes)) {
    rnk <- -1*log10(generes[[i]]$p_val_adj) * sign(generes[[i]]$avg_logFC)
    # compute the padj and fold change into the rank
    names(rnk) <- rownames(generes[[i]])
    scmappr[names(rnk),i] <- unname(rnk)
  }
  wilcoxon_rank_mat_t <- scmappr
  wilcoxon_rank_mat_t[is.infinite(wilcoxon_rank_mat_t) & wilcoxon_rank_mat_t < 0] <- min(wilcoxon_rank_mat_t[is.finite(wilcoxon_rank_mat_t)])
  wilcoxon_rank_mat_t[is.infinite(wilcoxon_rank_mat_t) & wilcoxon_rank_mat_t > 0] <- max(wilcoxon_rank_mat_t[is.finite(wilcoxon_rank_mat_t)])
  
  
  # generating the signature matrix developed by the rank of fold changes
  scmappr <- matrix(-300, nrow = length(genes_uni), ncol = length(names(generes)))
  rownames(scmappr) <- genes_uni
  colnames(scmappr) <- names(generes)
  for(i in 1:length(generes)) {
    rnk <- generes[[i]]$avg_logFC
    names(rnk) <- rownames(generes[[i]])
    scmappr[names(rnk),i] <- unname(rnk)
  }
  wilcoxon_rank_mat_or <- scmappr
  wilcoxon_rank_mat_or[is.infinite(wilcoxon_rank_mat_or) & wilcoxon_rank_mat_or < 0] <- min(wilcoxon_rank_mat_or[is.finite(wilcoxon_rank_mat_or)])
  wilcoxon_rank_mat_or[is.infinite(wilcoxon_rank_mat_or) & wilcoxon_rank_mat_or > 0] <- max(wilcoxon_rank_mat_or[is.finite(wilcoxon_rank_mat_or)])
  wilcoxon_rank_mat_or <- 2^wilcoxon_rank_mat_or # instead of the log2Fc as the output, ouput the odds ratio (only positive) that the gene is in that cell-type
  l <- list(pVal = wilcoxon_rank_mat_t, OR = wilcoxon_rank_mat_or, cellname = cell_name, topGenes = topGenes)
  # return the pVal signature matrix, odds ratio signautre matrix, cell names for that list, and topGenes
  return(l)
}

