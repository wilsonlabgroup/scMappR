#' Identify co-expressed cell-types.
#'
#' This function identifies genes with similar cell-type markers and if those markers are driving enrichment.
#'
#' This function takes significantly enriched cell-types from the single CT enrich before testing to see if the genes driving their enrichment are overlapping.
#' To save computational time and to not complete this with an incredible number of permutations, scMappR stops at overlapping 5 cell-types.
#'
#'
#' @rdname coEnrich
#' @name coEnrich
#'
#' @param sig A The number of combinations of significant cell-types to enrich.
#' @param gene_list_heatmap Signature matrix of inputted genes in heatmap and the cell-type preferences -- output of heatmap generation.
#' @param background_heatmap  Signature matrix of background matrix in heatmap and cell-type preferences -- output of heatmap generation.
#' @param study_name Name of the outputted table.
#' @param output_directory Name of the directory this table will be printed in.
#' 
#'
#' @return \code{coEnrich} Enrichment of cell-types that are expressed by the same genes, from 2-5 sets of cell-types. \cr
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
#'
#' @examples 
#' 
#' 
#' # load in signature matrices
#'  load("data/Preoptic_region_example.rda")
#' # data(Preoptic_region_example)
#'  Signature <- POA_Rank_signature
#'  genes <- rownames(Signature)[1:200]
#'  heatmap_test <- tissue_scMappR_custom( genes, signature_matrix = Signature,output_directory =  "scMappR_test")
#'  
NULL

#' @rdname coEnrich
#' @export
#' 
coEnrich <- function(sig, gene_list_heatmap, background_heatmap, study_name, outDir = output_directory) {
  # Internal
  # this function takes significantly enriched cell-types from the single CT enrich before testing to see
  # if the genes driving their enrichment are overlapping
  # to save computational time and to not complete this with an incredible number of permutations, scMappR stops at overlapping 5 cell-types
  # Args:
  # sig = preferences of significant functions
  # gene_list_heatmap = signature matrix (p-value) of your subsetted gene list
  # background_heatmap = signature matrix (p-value) of entire tissue
  # study_name = the name of the co-enrich table to be outputted
  # output_directory = the name of the directory where co-enrich will be completed
  # Returns:
  # Enrichment of cell-types that are expressed by the same genes, from 2-5 sets of cell-types.
  if(nrow(sig) > 5) {
    sig <- sig[1:5,]
  }
  l <- nrow(sig)
  multi_comps <- c()
  for(y in 2:l) {
    # for combinations of 2-# enriched cell types (max= 5)
    if(y < l) {
      comps <- utils::combn(sig$cell_type, y)
      # find combinations of cell-types
      co_up <- function(x) return(length(x[x>=1])==y)
      #score genes significantly enriched in both cell types = 1
      print(y)
      for(j in 1:ncol(comps)) {
        # for that comparison
        thecomps <- toupper(comps[,j]) # take comparison
        geneList_comb <- gene_list_heatmap$geneHeat # get gene list signature matrix
        colnames(geneList_comb) <- toupper(colnames(geneList_comb))
        geneList_comb1 <- geneList_comb[,which(colnames(geneList_comb) %in% thecomps)] # extract cell-types to be cenriched
        background_comb <- background_heatmap$geneHeat # take tissue signature matrix
        colnames(background_comb) <- toupper(colnames(background_comb))
        
        inter <- intersect(rownames(background_heatmap), rownames(geneList_comb))
        
        background_comb1 <- background_comb[!(rownames(background_comb) %in% inter), which(colnames(background_comb) %in% thecomps)] # remove enriched CT
  
        bin_aa <- apply(geneList_comb1,1,co_up) # see who is co-enriched in gene list
        name_in <- names(bin_aa)[bin_aa == T]
        aa <- sum(bin_aa) # number of co-enriched genes in list
        ab <- nrow(geneList_comb1) - aa # number non-co-enriched genes not in list
        ba <- sum(apply(background_comb1,1,co_up)) # number of co-enriched not in list
        bb <- nrow(background_comb1) - ba # number of non-coenriched genes not in list
        m <- matrix(c(aa, ba, ab, bb), nrow = 2) # fishers exact test
        fisherTest <- stats::fisher.test(m)
        OR <- fisherTest$estimate
        p <- fisherTest$p.value
        name <- paste0(thecomps,collapse=":" )
        row <- c(name,p,OR, paste0(name_in, collapse = ","))
        multi_comps <- rbind(multi_comps, row)
      }
    }
    if(y == l ) { # if you're on the the max nmber of co-enriched celltypes you do the same thing but without the combination step since it's the max number
      co_up <- function(x) return(length(x[x>=1])==y)
      
      thecomps <- sig$cell_type
      geneList_comb <- gene_list_heatmap$geneHeat
      colnames(geneList_comb) <- toupper(colnames(geneList_comb))
      geneList_comb1 <- geneList_comb[,thecomps]
      background_comb <- background_heatmap$geneHeat
      colnames(background_comb) <- toupper(colnames(background_comb))
      
      inter <- intersect(rownames(background_heatmap), rownames(geneList_comb))
      background_comb1 <- background_comb[!(rownames(background_comb) %in% inter), which(colnames(background_comb) %in% thecomps)] # remove enriched CT
      
      bin_aa <- apply(geneList_comb1,1,co_up)
      name_in <- names(bin_aa)[bin_aa == T]
      aa <- sum(bin_aa)
      ab <- nrow(geneList_comb1) - aa
      ba <- sum(apply(background_comb1,1,co_up))
      bb <- nrow(background_comb1) - ba
      m <- matrix(c(aa, ba, ab, bb), nrow = 2)
      fisherTest <- stats::fisher.test(m)
      OR <- fisherTest$estimate
      p <- fisherTest$p.value
      name <- paste0(thecomps,collapse=":" )
      row <- c(name,p,OR, paste0(name_in, collapse = ","))
      multi_comps <- rbind(multi_comps, row)
    }
  }
  colnames(multi_comps) <- c("cell_types", "p_val", "OR", "genes")
  write.table(multi_comps, file = paste0(outDir, "/",study_name, "cell_co_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
  multi_comps <- read.table(file = paste0(outDir, "/",study_name, "cell_co_preferences.tsv"), as.is=T,header=T, sep = "\t")
  multi_comps$pFDR <- p.adjust(multi_comps$p_val, "fdr")
  write.table(multi_comps, file = paste0(outDir, "/",study_name, "cell_co_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
  return(multi_comps)
}

