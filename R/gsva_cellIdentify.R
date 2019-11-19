#' Cell-type naming with GSVA
#'
#' This function computes the mean expression of every cell-typ before predicting the most likely cell-type using the GSVA method.s
#'
#' This function inputs a Seurat object and uses the average normalized expression of each gene in each cluster to identify cell-types using the gsva method.
#'
#' @rdname gsva_cellIdentify
#' @name gsva_cellIdentify
#'
#' @param pbmc Processed seurat object without named cells
#' @param theSpecies "human" or "mouse" -- it will determine which CT marker database to use -- there are some differences
#' @param naming_preference Once top CT markers are identified, naming_preferences will then extract CT markers within a more appropriate tissue type etc.
#' @param rda_path Path to precomputed cell-type gmt files (rda objects).
#' 
#'
#' @return \code{gsva_cellIdentify} A list containing the top cell-type marker for a cell-type using the panglao dataset as well as the cellMarker dataset. \cr
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
#' load("~/scMappR/data/cell_process_example.rda")
#' #data(cell_process_example)
#' toProcess <- list(example = sm)
#' tst1 <- process_from_count(toProcess, "testProcess")
#' cellnames <- gsva_cellIdentify(tst1, "mouse", "brain", "~/scMappR/data")
#' 
#'  
#' @export
#' 
gsva_cellIdentify <- function(pbmc, theSpecies, naming_preference, rda_path = "") {
  # this function inputs a Seurat object and uses the average normalized expression of each gene
  # in each cluster to identify cell-types using the gsva method
  # Args: 
  #pbmc: processed seurat object
  # theSpecies: "human" or "mouse" -- it will determine which CT marker database to use -- there are some differences
  # naming_preference: once top CT markers are identified, naming_preferences will then extract CT markers within a more appropriate tissue type etc.
  # rda_path Path to precomputed cell-type gmt files (rda objects).
  # Returns: 
  # A list containing the top cell-type marker for a cell-type using the panglao dataset as well as the cellMarker dataset
  
  avg_expr <- Seurat::AverageExpression(pbmc)
  # identify average expression of clusters  
  # panglao
  print(theSpecies)
  thefiles <- list.files(path = rda_path, "_cell_markers.rda")

  if(theSpecies == "human") { # load cell marker database
    if(length(thefiles) == 0) {
      warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
      #
      datafile <- "human_cell_markers.rda"
      metafile <- paste0(datafile)
      url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      #
    } else {
      load(paste0(rda_path,"/human_cell_markers.rda"))
    }
  } else {
    if(length(thefiles) == 0) {
      warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
      #
      datafile <- "mouse_cell_markers.rda"
      metafile <- paste0(datafile)
      url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      #
    } else{
    load(paste0(rda_path,"/mouse_cell_markers.rda"))
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
    top5_pang[[i]] <-  round(sort(gbm_pang[,i], decreasing = T)[1:5],2)
    top5_cm[[i]] <- round(sort(gbm_cellmarker[,i], decreasing = T)[1:5],2)
  }
  
  if(length(thefiles) == 0) {
    warning(paste0("Cell-marker preferences are not present in ", rda_path, " downloading and loading data."))
    #
    datafile <- "cell_preferences_categorized.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
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
  
  l <- list(cellMarker = cm_top, panglao = pang_top)
  return(l)
}

