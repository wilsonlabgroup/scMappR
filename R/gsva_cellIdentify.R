#' Cell-type naming with GSVA
#'
#' This function computes the mean expression of every cell-type before predicting the most likely cell-type using the GSVA method.
#'
#' This function inputs a Seurat object and uses the average normalized expression of each gene in each cluster to identify cell-types using the GSVA method.
#'
#' @rdname gsva_cellIdentify
#' @name gsva_cellIdentify
#'
#' @param pbmc Processed seurat object without named cells.
#' @param theSpecies "human" or "mouse" -- it will determine which CT marker database to use -- there are some differences.
#' @param naming_preference Once top CT markers are identified, naming_preferences will then extract CT markers within a more appropriate tissue type.
#' @param rda_path Path to pre-computed cell-type .gmt files (rda objects).
#' @param toSave If scMappR is allowed to write files and directories.
#'
#' @return List with the following elements:
#' \item{cellMarker}{Most likely cell-types predicted from cellMarker database.}
#' \item{panglao}{Most likely cell-types predicted from panglao database.}
#' \item{avg_expression}{Average expression of each gene in each cell-type.}
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
#' tst1 <- process_from_count(toProcess, "testProcess")
#' cellnames <- gsva_cellIdentify(tst1, theSpecies = "mouse",
#'  naming_preference = "brain", rda_path = "")
#' 
#' }
#'  
#' @export
#' 
gsva_cellIdentify <- function(pbmc, theSpecies, naming_preference = -9, rda_path = "", toSave = FALSE) {
  # this function inputs a Seurat object and uses the average normalized expression of each gene
  # in each cluster to identify cell-types using the gsva method
  # Args: 
  #pbmc: processed seurat object
  # theSpecies: "human" or "mouse" -- it will determine which CT marker database to use -- there are some differences
  # naming_preference: once top CT markers are identified, naming_preferences will then extract CT markers within a more appropriate tissue type etc.
  # rda_path Path to precomputed cell-type gmt files (rda objects).
  # Returns: 
  # A list containing the top cell-type marker for a cell-type using the panglao dataset as well as the cellMarker dataset
  not_seurat <- class(pbmc)[1] == "Seurat"
  if((not_seurat[1] == FALSE)[1]) {
    stop("pbmc must be of class 'Seurat'")
  }

  naming_preferences <- c("brain", "epithelial", "endothelial", "blood", "connective","eye", "epidermis", "Digestive", "Immune", "pancreas", "liver", "reproductive", "kidney", "respiratory") 
  if(!naming_preference %in% naming_preferences) {
    if(naming_preference != -9) {
      message("Naming preference options")
    message(naming_preferences)
    stop("Naming preferences not in options (case sensitive) and isn't a non-choice (-9), please try again.")
    }
  }
  if(!(theSpecies %in% c("human", "mouse"))) {
    if(theSpecies != -9) {
      stop("theSpecies is not 'human' 'mouse' or '-9' (case sensitive), please try again with this filled.")
    }
  }
  
  if(!is.character(rda_path)) {
    stop("rda_path must be of class character.")
  }
  
  if(!is.logical(toSave)) {
    stop("toSave must be of class logical.")
  }
  
  
  avg_expr <- Seurat::AverageExpression(pbmc)
  # identify average expression of clusters  
  # panglao
  message(theSpecies)
  thefiles <- list.files(path = rda_path, "_cell_markers.rda")
  gmt_list <- "" # to help with the CRAN warning
  cell_preference_final <- ""
  if(theSpecies == "human") { # load cell marker database
    if(length(thefiles) == 0) {
      warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
      #
      datafile <- "human_cell_markers.rda"
      metafile <- paste0(datafile)
      url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      #
      gmt_both <- gmt_list$gmt_both
      gmt_cellmarker <- gmt_list$gmt_cellmarker
      gmt_gobp <- gmt_list$gmt_gobp
      gmt_panglao <- gmt_list$gmt_panglao
      gmt_subtype <- gmt_list$gmt_subtype
    } else {
      load(paste0(rda_path,"/human_cell_markers.rda"))
      gmt_both <- gmt_list$gmt_both
      gmt_cellmarker <- gmt_list$gmt_cellmarker
      gmt_gobp <- gmt_list$gmt_gobp
      gmt_panglao <- gmt_list$gmt_panglao
      gmt_subtype <- gmt_list$gmt_subtype
    }
  } else {
    if(length(thefiles) == 0) {
      warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
      #
      datafile <- "mouse_cell_markers.rda"
      metafile <- paste0(datafile)
      url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      gmt_both <- gmt_list$gmt_both
      gmt_cellmarker <- gmt_list$gmt_cellmarker
      gmt_gobp <- gmt_list$gmt_gobp
      gmt_panglao <- gmt_list$gmt_panglao
      gmt_subtype <- gmt_list$gmt_subtype
      #
    } else{
    load(paste0(rda_path,"/mouse_cell_markers.rda"))
      gmt_both <- gmt_list$gmt_both
      gmt_cellmarker <- gmt_list$gmt_cellmarker
      gmt_gobp <- gmt_list$gmt_gobp
      gmt_panglao <- gmt_list$gmt_panglao
      gmt_subtype <- gmt_list$gmt_subtype
    }
  }
  gmt <- gmt_both
  gbm <- GSVA::gsva(as.matrix(avg_expr$RNA), gmt, mx.diff=FALSE, verbose=FALSE, parallel.sz=1, min.sz= 5)
  # complete GSVA of the average expression of each cell-types
  gbm_pang <- gbm[grep("panglao", rownames(gbm)),] # cell-types from panglao
  gbm_cellmarker <- gbm[-grep("panglao", rownames(gbm)),] # cell-types from cellmarker
  
  # get the top 5 cell-types from each dataset
  top5_pang <- list()
  top5_cm <- list()
  for(i in 1:ncol(gbm)) {
    top5_pang[[i]] <-  round(sort(gbm_pang[,i], decreasing = TRUE)[1:5],2)
    top5_cm[[i]] <- round(sort(gbm_cellmarker[,i], decreasing = TRUE)[1:5],2)
  }
  
  if(length(thefiles) == 0) {
    warning(paste0("Cell-marker preferences are not present in ", rda_path, " downloading and loading data."))
    #
    datafile <- "cell_preferences_categorized.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
  } else {
    
    load(paste0(rda_path,"/cell_preferences_categorized.rda"))
  }
  
  
  top_ct <- function(x) { 
    # If you have a cell-type naming preference and cell-types that lie within that preference, those are extracted
    # the top cell-type for 
    if(naming_preference != -9) {
      
      # if they have the prefered tissue or cell-type to pick from
      # prioritize those
      mypref <- cell_preference_final[[naming_preference]]
      top5_ct_pref <- names(x)
      outcell_preference <- which(top5_ct_pref %in% mypref)
      if(length(outcell_preference) > 0 ) {
        y <- x[outcell_preference[1]]
        return(paste0(names(y[1]),"_",unname(y[1])))
      } else {
        return(paste0(names(x[1]),"_",unname(x[1])))
      }
      
    } else {
      return(paste0(names(x[1]),"_",unname(x[1])))
    }
  }
  # return the top cell-type from each dataset in a manner is that cognisant of cell-type preferences
  cm_top <- lapply(top5_cm, top_ct)
  pang_top <- lapply(top5_pang, top_ct)
  ave_expr_RNA <- avg_expr$RNA
  l <- list(cellMarker = cm_top, panglao = pang_top, avg_expression = ave_expr_RNA)
  return(l)
}

