#' Fisher's Exact Cell-Type Identification.
#'
#' This function uses the cellmarker and panglao datasets to identify cell-type DEGs
#' 
#' Complete a fishers exact test of an input list of genes against a gene set saved in an *.RData object.
#' The RData object is storing a named list of genes called "gmt"
#'
#' @rdname cellmarker_enrich
#' @name cellmarker_enrich
#'
#' @param gene_list: a character vector of gene symbols with the same designation (e.g. mouse symbol - mouse, human symbol - human) as the gene set database
#' @param p_thresh: the Fisher's test cutoff for a cell-marker to be enriched.
#' @param gmt: either a path to an rda file containing an object called "gmt", which is a named list where each element of the list is a vector of gene symbols website for more detail on the file type (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats).
#' @param the gmt list may also be inputted.  
#' @param fixed_length Estimated number of genes in your background.
#' @param min_genes Minimum number of genes in the cell-type markers.
#' @param max_genes Maximum number of genes in a cell-type marker.
#' @param isect_size Number of genes in your list and the cell-type.
#'
#' @return \code{cellmarker_enrich} Gene set enrichment of cell-types on your inputted gene list. \cr
#'
#'
#' @examples 
#' \notrun {
#' 
#' # load in signature matrices
#' load("data/Preoptic_region_example.rda")
#' # data(Preoptic_region_example)
#' Signature <- POA_Rank_signature
#'  rowname <- get_gene_symbol(Signature)
#'  rownames(Signature) <- rowname$rowname
#'  genes <- rownames(Signature)[1:100]
#'  
#'  # Assuming mouse_cell_markers.rda is in you "~/Documents/scMappR/data" directory
#'  gmt1 <- "~/Documents/scMappR/data/mouse_cell_markers.rda"
#'  gmt <- gmt_panglao
#'  enriched <- cellmarker_enrich(genes, 0.05, gmt = gmt)
#'  }
NULL

#' @rdname cellmarker_enrich
#' @export
#' 
cellmarker_enrich <- function(gene_list, p_thresh, gmt = "cellmarker_list.Rdata", fixed_length = 13000, min_genes = 5, max_genes = 3000, isect_size = 3) {
  # complete a fishers exact test of an input list of genes against a gene set saved in an *.RData object
  # The RData object is storing a named list of genes called "gmt"
  
  # Args:
  # gene_list: a list of gene symbols with the same designation (e.g. mouse symbol - mouse, human symbol - human) as the gene set database
  # p_thresh: the Fisher's test cutoff for a cell-marker to be enriched
  # gmt: either a path to an rda file containing an object called "gmt", which is a named list where each element of the list is a vector of gene symbols
  # website for more detail on the file type (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)
  # the gmt list may also be inputted.  
  # fixed_length = estimated number of genes in your background
  # min_genes: minimum number of genes in the cell-type markers
  # max_genes: maximum number of genes in a cell-type marker
  # isect_size: number of genes in your list and the cell-type
  
  # Returns: 
  # gene set enrichment of cell-types on your inputted gene list
  
  
  toTest_n <- gene_list
  if(class(gmt) == "character") {
    load(gmt) # load the gmt object
  }
  theRows <- c()
  for(i in 1:length(gmt)) {
    
    test <-gmt[[i]]
    if(length(test) < min_genes | length(test) > max_genes) {
      next
    }
    ints <- intersect(toTest_n, test) # intersect your inputted genes with the study
    tbt = matrix(c(length(ints), length(toTest_n)-length(ints),  length(test) ,fixed_length-length(test)), nrow=2)
    # make a 2x2 fisher's matrix
    rownames(tbt) <- c("Exposure_yes", "Exposure_no")
    colnames(tbt) <- c("outcome_yes", "outcome_no")
    f_out <- fisher.test(tbt, B = 1e9, alternative = "greater") # fishers exact test for over-representation
    p <- f_out$p.value
    OR <- f_out$estimate[[1]]
    if(p < 0.05 & OR < 1) {
      p = 1 # assume overrepressentation
    }
    term_size <- length(test)
    intersect_size <- length(ints)
    input_length <- length(toTest_n)
    proportion_success <- signif(intersect_size / term_size,2)
    interGenes <- paste0(ints, collapse = ",")
    theName <- names(gmt)[i]
    theRow <- c(theName, p, term_size, intersect_size, input_length, interGenes) 
    theRows <- rbind(theRows, theRow)
  }
  colnames(theRows) <- c("name","p", "term_size", "intersect_size", "input_length", "genes") # save study as gene set enrichment output
  df_theRows <- as.data.frame(theRows)
  # convert cell-type enrichment matrix 
  df_theRows$p <- toNum(df_theRows$p)
  fdr <- p.adjust(df_theRows$p, "fdr")
  bonf <- p.adjust(df_theRows$p, "bonferroni")
  df_theRows$fdr <- fdr
  df_theRows$bonf <- bonf
  df_theRows$intersect_size <- toNum(df_theRows$intersect_size)
  df_theRows_order <- df_theRows[order(df_theRows$p),]
  
  
  df_theRows_order <- df_theRows_order[df_theRows_order$intersect_size > isect_size-1,]
  
  sig_studies <- df_theRows_order[df_theRows_order$p < p_thresh,]
  return(sig_studies)
}

