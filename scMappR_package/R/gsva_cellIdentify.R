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
#'
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
NULL

#' @rdname gsva_cellIdentify
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
  library(GSVA)
  avg_expr <- AverageExpression(pbmc)
  # identify average expression of clusters  
  # panglao
  print(theSpecies)
  if(theSpecies == "human") { # load cell marker database
    load(paste0(rda_path,"/human_cell_markers.rda"))
  } else {
    load(paste0(rda_path,"/mouse_cell_markers.rda"))
  }
  gmt <- gmt_both
  gbm <- gsva(as.matrix(avg_expr$RNA), gmt, mx.diff=FALSE, verbose=FALSE, parallel.sz=1, min.sz= 5)
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
  
  top_ct <- function(x) { 
    # If you have a cell-type naming preference and cell-types that lie within that preference, those are extracted
    # the top cell-type for 
    if(naming_preference != -9) {
      load(paste0(rda_path,"/cell_preferences_categorized.rda"))
      #data(cell_preferences_categorized)
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

