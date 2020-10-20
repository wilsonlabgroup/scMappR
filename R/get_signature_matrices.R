#' Get signature matrices.
#'
#' This function downloads and returns signature matrices and associated cel-type labels from the scMappR_data repo.
#'
#'
#' @rdname get_signature_matrices
#' @name get_signature_matrices
#'
#' @param type a character vector that can be 'all', 'pVal', or 'OR'
#'
#' @return \code{get_signature_matrices} Returns the signature matrices currently stored in scMappR_Data. Associated cell-type labels from different methods for each signature matrix is also provided.\cr
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases ks.test
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
#' signatures <- get_signature_matrices(type = "all")
#'  
#' }
#'  
#' @export

get_signature_matrices <- function(type = "all") {

  type_check <- !all(is.character(type),length(type) == 1, type %in% c("all", "pVal", "OR", "label"))
  if(type_check[1]) {
    stop("the 'type' variable must be a character, specifically 'all', 'pVal', 'OR', or 'label'")
  }
  
  RankValueSignature <- "" # empty for cran
  OddsRatioSignature <- "" # 
  celltypeLabels <- "" # 
  
  if(type == "all") {
    
    ## The signature matrix (p-value)
    datafile <- "Signature_matrices_Pval.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    pVal <- RankValueSignature
    
    ## The signature matrix (odds-ratio)
    datafile <- "Signature_matrices_OR.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    OR <- OddsRatioSignature
    
    ## The matrix of cell-type labels
    datafile <- "Signature_matrices_cellLabel.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    labels <- celltypeLabels
    
    l <- list(pVal = pVal, OR = OR, labels = labels )
    return(l)
  }
  if(type == "pVal") {
    ## The signature matrix (p-value)
    datafile <- "Signature_matrices_Pval.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    pVal <- RankValueSignature
   return(pVal)
  }
  
  if(type == "OR") {
    ## The signature matrix (p-value)
    datafile <- "Signature_matrices_Pval.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    OR <- OddsRatioSignature
    return(OR)
  }
  if(type == "label") {
    ## The matrix of cell-type labels
    datafile <- "Signature_matrices_cellLabel.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    labels <- celltypeLabels
    return(labels)
  }
  
}
