#' Gene List Visualization and Enrichment with Custom Signature Matrix
#'
#' This function visualizes signature matrix, clusters subsetted genes, completes enrichment of individual cell-types and co-enrichment.
#' 
#' This function is roughly the same as tissue_scMappR_internal, however now there is a custom signature matrix.
#' it generates a heatmap of the signature matrix and your inputted gene list as well as single cell-type and 
#' co-celltype enrichment.
#'
#'
#' @rdname tissue_scMappR_custom
#' @name tissue_scMappR_custom
#'
#' @param gene_list A list of gene symbols matching that of the signature_matrix -- any gene symbol is acceptable.
#' @param signature_matrix Precomputed signature matrix with matching gene names.
#' @param output_directory Directory made containing output of functions.
#' @param gene_cutoff Value cutoff (generally rank := log10(Padj)) for a gene to be considered a marker.
#' @param is_pvalue If signature matrix is p-value before rank is applied (not recommended ) (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#'
#' @return \code{tissue_scMappR_custom} A list containing the entire signature matrix, the matrix subsetted for your genes, enrichment of each cell-type, and co-enrichment. \cr
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
#' @importFrom gProfileR gprofiler
#' @importFrom gprofiler2 gost
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#'
#' @examples 
#' 
#' 
#' # load in signature matrices
#' data(POA_example)
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' Signature <- POA_Rank_signature
#' genes <- rownames(Signature)[1:60]
#' heatmap_test <- tissue_scMappR_custom( genes, signature_matrix = Signature,
#'                                       output_directory =  "scMappR_test", toSave = FALSE)
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
  
  if(class(gene_list) != "character") {
    stop("gene_list must be of class character.")
  }
  if(class(signature_matrix) != "data.frame" & class(signature_matrix) != "matrix") {
    stop("signature_matrix must be of class data.frame or matrix.")
  }
  if(class(output_directory) != "character") {
    stop("output_directory must be of class character.")
  }
  if(class(gene_cutoff) != "numeric") {
    stop("gene_cutoff must be of class numeric." )
  }
  if(!(any(is.logical(toSave), is.logical(is_pvalue)))) {
    stop("toSave and is_pvalue must be of class logical.")
  }
  
  
  outDir <- paste0(output_directory)
  study_names <- outDir
  if(toSave == TRUE) {
    dir.create(outDir)
  } else {
    warning("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.")
    print("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.", quote = FALSE)
    
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
    utils::write.table(singleCTpreferences, file = paste0(outDir,"/",outDir, "_celltype_preferences.tsv"), quote= FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
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

