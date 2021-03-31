#' Identify all cell-type markers 
#'
#' Takes processed Seurat matrix and identifies cell-type markers with FindMarkers in Seurat.
#'
#' Internal: This function runs the FindMarkers function from Seurat in a loop, will use the Seurat v2 or Seurat v3 object after identifying which Seurat object is inputted. 
#' It then takes the output of the FindMarkers and puts it in a list, returning it.
#' 
#' @rdname seurat_to_generes
#' @name seurat_to_generes
#'
#' @param pbmc Processed Seurat object.
#' @param test statistical test for calling CT markers -- must be in Seurat.
#'  
#'
#' @return \code{seurat_to_generes} A list of genes where their over-representation in the i'th cell-type is computed. Each element contains the gene name, adjusted p-value, and the log2Fold-Change of each gene being present in that cell-type. \cr
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text geom_bar theme_classic xlab ylab scale_fill_manual element_line
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases ks.test dist shapiro.test mad
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#' @importFrom ADAPTS estCellPercent
#' @importFrom reshape melt
#'
#' @examples
#' \donttest{  
#'  
#' data(sm)
#' toProcess <- list(example = sm)
#' tst1 <- process_from_count(toProcess, "testProcess", theSpecies = "mouse")
#' generes <- seurat_to_generes(tst1)
#' }
#' 
#' @export



seurat_to_generes <- function(pbmc, test = "wilcox"){
  # Internal: This function runs the FindMarkers function from seurat in a loop, will use the seurat v2 or seurat v3 object 
  # after identifying which seurat object is inputted. It then takes the output of the FindMarkers and puts it in a list, returning it
  # Args:
  # pbmc -- a processed suerat object
  # Returns:
  # A list of genes where their over-representation in the i'th cell-type is computed. Each element contains the gene name, adjusted p-value, and the log2FC of each gene being present in that cell-type.
  isSeurat <- class(pbmc)[1] == "Seurat"
  if(isSeurat[1] == FALSE) {
    stop("pbmc object must be of class Seurat.")
  }
  
  id <- try(pbmc@ident, silent = TRUE)
  # try for seurat v2
  
  generes <- list()
  count <-1
  isError <- class(id)[1] == "try-error"
  if(isError == FALSE) {
    # if it's V2
    for(i in sort(unique(pbmc@ident))) {
      
      de_genes <- try(Seurat::FindMarkers(pbmc, ident.1 = i, test.use = test))
      worked <- any(class(de_genes)[1] == "try-error", is.null(de_genes))
      if(worked[1]) {
        
        next
      }
      
      generes[[count]] <- de_genes
      names(generes)[count] <- i
      count <- count+1
      message("DE")
      message(i)
    }
    
    
  }
  for(i in sort(unique(pbmc@active.ident))) {
    # if object is Seurat V3
    de_genes <- try(Seurat::FindMarkers(pbmc, ident.1 = i, test.use = test))
    worked <- any(class(de_genes)[1] == "try-error", is.null(de_genes))
    
    if(worked[1]) {
      
      next
    }
    
    generes[[count]] <- de_genes
    names(generes)[count] <- i
    count <- count+1
    message("DE")
    message(i)
  }
  if(is.null(generes[[1]]$avg_logFC)) {
    message("Seurat V4 or later was used to identify cell-type markers, adding 'avg_logFC' column. It has the same data as avg_log2FC but is compatible with downstream functions in scMappR.")
    for(z in 1:length(generes))
      generes[[z]]$avg_logFC <- generes[[z]]$avg_log2FC
  }
  
  return(generes)
}

