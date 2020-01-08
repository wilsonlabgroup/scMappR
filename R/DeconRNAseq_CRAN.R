##' DeconRNASeq CRAN compatible
#' 
#' This function runs DeconRNAseq with default parameters such that it is compatible with CRAN and scMappR
#'
#' This is the exact same function as the primary function in the Bioconductor package, DeconRNAseq (PMID: 23428642)
#' except it is now compatible with CRAN packages. 
#'
#' @rdname DeconRNAseq_CRAN
#' @name DeconRNAseq_CRAN
#'
#' @param datasets Normalized RNA-seq dataset
#' @param signatures Signature matrix of odds ratios
#' @param proportions If cell-type proportion is already inputted - always NULL for scMappR
#' @param checksig Check to see if plotting is significant - always false for scMappR
#' @param known.prop If proportions were known - always false for scMappR
#' @param use.scale Scale and center value - always TRUE for scMappR
#' @param fig Make figures - always FALSE for scMappR
#' 
#' @return \code{DeconRNAseq_CRAN} Estimated cell-type proportions with DeconRNAse \cr
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
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#'
#' @examples 
#' 
#' data(PBMC_example)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' out <- DeconRNAseq_CRAN(as.data.frame(bulk_normalized), as.data.frame(odds_ratio_in))
#'                                      
#' @export
#' 

DeconRNAseq_CRAN <- function (datasets, signatures, proportions = NULL, checksig = FALSE,  known.prop = FALSE, use.scale = TRUE, fig = FALSE) {
  # datasets normalized RNA-seq dataset
  # Signatures signature matrix of odds ratios
  # proportions if cell-type proportion is already inputted - always NULL for scMappR
  # checksig check to see if plotting is significant - always false for scMappR
  # known.prop if proportions were known - always false for scMappR
  # use.scale Scale and center value - always TRUE for scMappR
  # fig Make figures - always FALSE for scMappR
  if (is.null(datasets)) 
    stop(" Missing the mixture dataset, please provide a tab-delimited text file for mixture samples.")
  if (is.null(signatures)) 
    stop(" Missing the signature dataset, please provide a tab-delimited text file for pure tissue/cell types.")
  if (is.null(proportions) && known.prop) 
    stop(" Missing the known proprotions, please provide a tab-delimited text file containing known fractions for pure tissue/cell types.")
  x.signature <- signatures
  x.data <- datasets
  if (is.data.frame(x.signature) == FALSE) 
    stop("signature datasets must be a dataframe")
  if (sum(is.na(x.signature)) > 0) 
    stop("signature data cannot have NAs. please exclude or impute missing values.")
  if (is.data.frame(x.data) == FALSE) 
    stop("mixture datasets must be a dataframe")
  if (sum(is.na(x.data)) > 0) 
    stop("mixture data cannot have NAs. please exclude or impute missing values.")
  numofg <- nrow(x.signature)
  Numofx <- ncol(x.signature)
  if (numofg < Numofx) 
    stop("The number of genes is less than the number of cell types, which means less independent equations than unknowns.")
  x.data.temp <- pcaMethods::prep(x.data, scale = "none", center = TRUE)
  x.data.pca <- pcaMethods::pca(x.data.temp, method = "svd", center = FALSE, 
                                nPcs = Numofx)
  out.pca <- ""
  Var <- pcaMethods::R2cum(x.data.pca)
  numofmix <- order(Var > 0.99, decreasing = T)[1]
  if (checksig && numofg >= 40) {
    step <- seq(20, numofg, by = 20)
    sig.cond <- sapply(step, function(x) kappa(scale(x.signature[1:x, 
                                                                 ])))
  }
  common.signature <- rownames(x.signature) %in% rownames(x.data)
  common.data <- rownames(x.data) %in% rownames(x.signature)
  x.data <- x.data[common.data, ]
  x.signature <- x.signature[common.signature, ]
  x.subdata <- x.data[rownames(x.signature), ]
  Numofx <- ncol(x.signature)
  if (use.scale) {
    AA <- scale(x.signature)
  }
  else {
    AA <- x.signature
  }
  EE <- rep(1, Numofx)
  FF <- 1
  GG <- diag(nrow = Numofx)
  HH <- rep(0, Numofx)
  out.all <- c()
  for (i in colnames(x.subdata)) {
    BB <- x.subdata[, i]
    if (use.scale) {
      BB <- scale(BB)
    }
    out <- limSolve::lsei(AA, BB, EE, FF, GG, HH)
    out.all <- rbind(out.all, out$X)
  }
  mean.rmse <- 0
  rmse <- c()
  
  
  return(list(out.all = out.all, out.pca = out.pca))
  
}
