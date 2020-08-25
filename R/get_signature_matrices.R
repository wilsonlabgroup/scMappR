#' Get signature matrices.
#'
#' This function downloads and returns signature matrices from the scMappR_data repo.
#'
#'
#' @rdname get_signature_matrices
#' @name get_signature_matrices
#'
#' @param type a character vector that can be 'both', 'pVal', or 'OR'
#'
#' @return \code{get_signature_matrices} Returns a list of signature matrices currently stored in scMappR_Data. \cr
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#'
#' @examples
#' \donttest{
#'  
#' signatures <- get_signature_matrices(type = "both")
#'  
#' }
#'  
#' @export

get_signature_matrices <- function(type = "both") {

  type_check <- !all(is.character(type),length(type) == 1, type %in% c("both", "pVal", "OR"))
  if(type_check[1]) {
    stop("the 'type' variable must be a character, specifically 'both', 'pVal', or 'OR'")
  }
  
  RankValueSignature <- "" # empty for cran
  OddsRatioSignature <- "" # 
  
  if(type == "both") {
    
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
    
    l <- list(pVal = pVal, OR = OR)
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
  
}
