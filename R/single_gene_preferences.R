#'  Single cell-type gene preferences
#'
#'  Measure enrichment of individual cell-types in a signautre matrix
#'  
#' Internal function as part of tissue_scMappR_internal().
#' This function takes genes preferentially expressed within your gene list, each cell-type
#' and the background (i.e. all genes within the signature matrix) before completing identifying
#' the cell-type specific enrichment of the inputted gene list on each cell type.
#' This function then returns a table describing the cell-type enrichments (p-value and odds ratio)
#' of each cell-type
#'
#'
#' @rdname single_gene_preferences
#' @name single_gene_preferences
#'
#' @param hg_short A list with two objects: a "preferences" and a "genesIn". preferences is a list of gene symbols over-represented ine ach cell-type and genesIn were all the inputted genes.
#' @param hg_full The same as hg_short but fore very gene in the signature matrix.
#' @param study_name = name of output table
#' @param outputDir = directory where table is outputted
#' 
#'
#' @return \code{single_gene_preferences} A gene-set enrichment table of individual cell-type enrichment. \cr
#'
#'
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
#' @export
#' 
single_gene_preferences <- function(hg_short, hg_full, study_name, outDir = output_directory) {
  
  # Internal function as part of tissue_scMappR_internal()
  # This function takes genes preferentially expressed within your gene list, each cell-type
  # and the background (i.e. all genes within the signature matrix) before completing identifying
  # the cell-type specific enrichment of the inputted gene list on each cell type.
  # this function then returns a table describing the cell-type enrichments (p-value and odds ratio)
  # of each cell-type
  
  # Args:
  # hg_short = list of inputted genes, then inputted genes sorted to cell-type
  # hg_full = list of background genes, then background genes sorted to cell-type
  # study_name = name of output table
  # outputDir = directory where table is outputted
  # Returns:
  # A gene-set enrichment table of individual cell-type enrichment
  
  fP <- hg_full$preferences # your gene list sorted into each cell-type (if it's over-expressed)
  sP <- hg_short$preferences # all of the cell-type markers in the background
  in_length <- length(hg_short$genesIn) # the number of input genes
  inmarker_length <- length(hg_full$genesIn) # the number of background genes
  uni <- length(unique(c(hg_short$genesIn, hg_full$genesIn)))
  
  pref <- c()
  for(cell in 1:length(sP)) { # for every cell type
    theGenes <- unique(sP[[cell]]) # take unique markers
    theGenes_p <- paste0(theGenes,collapse = ",")
    aa <- length(unique(sP[[cell]])) # cell-type marker in gene list
    ab <- in_length - aa # gene list not cell-type marker
    ba <- length(unique(fP[[cell]])) # cell-type marker not gene lsit
    bb <-  uni - ba # background without that cell-type
    m <- matrix(c(aa, ba, ab, bb), nrow = 2) # fisher's test
    ftest <- fisher.test(m)
    p <- ftest$p.value
    OR <- ftest$estimate
    nSP <- names(sP)[cell]
    P <- c(nSP, p, OR, theGenes_p)
    pref <- rbind(pref, P) # combine the summary statistics
  }
  colnames(pref) <- c("cell_type", "p_val", "Odds_Ratio", "genes")
  # build the table and p-adjust 
  pref <- as.data.frame(pref)
  write.table(pref, file = paste0(outDir, "/",study_name, "cell_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
  pref <- read.table(file = paste0(outDir, "/",study_name, "cell_preferences.tsv"), as.is=T,header=T, sep = "\t")
  pref$pFDR <- p.adjust(pref$p_val, "fdr")
  write.table(pref, file = paste0(outDir, "/",study_name, "cell_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
  
  return(pref)
}

