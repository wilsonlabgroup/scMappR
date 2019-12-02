#' Consensus cell-type naming (Fisher's Exact)
#'
#' This function completes the Fisher's exact test cell-type naming for all cell-types
#' 
#' Fisher's exact test method of cell-type identification using the Panglao and CellMarker databases. It extracts significant pathways (pFDR < 0.05).
#' Then, if naming_preference != -9, it will extract the enriched cell-types within the cell-types identified with the naming preferences option.
#' Generally, this method seems to be biased to cell-types with a greater number of markers.
#'
#' @rdname human_mouse_ct_marker_enrich
#' @name human_mouse_ct_marker_enrich
#'
#' @param gene_lists A named list of vectors containing cell-type markers (mouse or human gene-symbols) which will be named as a cell-type via the fish'er sexact test method.
#' @param theSpecies The species of the gene symbols: human or mouse.
#' @param cell_marker_path If local, path to Cell-Type marker rda files, otherwise, we will try to download datafiles.
#' @param naming_preference Either -9 if there is no expected cell-type or one of the categories from get_naming_preference_options(). This is useful if you previously have an idea of which cell-type you were going to enrich.
#'
#' @return \code{human_mouse_ct_marker_enrich} list: The top cell-type for each cell-type markers list, as well as a matrix the top 5 most likely cell-type makers for each gene list (with OR and P-val) \cr
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
#' @import grDevices
#' @import gProfileR
#'
#' @examples 
#' \donttest{
#' 
# load in signature matrices
#' data(Preoptic_Area)
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' Signature <- POA_Rank_signature
#' rowname <- get_gene_symbol(Signature)
#' rownames(Signature) <- rowname$rowname
#' genes <- rownames(Signature)[1:100]
#' enriched <- human_mouse_ct_marker_enrich(gene_lists = genes, theSpecies = "mouse", 
#'                                          cell_marker_path = "", naming_preference = "brain")
#'  }
#' @export
#' 
human_mouse_ct_marker_enrich <- function(gene_lists, theSpecies = "human",cell_marker_path = "~/scMappR/data", naming_preference = -9) {
  # Fisher's exact test method of cell-type identification using the Panglao and CellMarker databases. It extracts significant pathways (pFDR < 0.05).
  # Then, if naming_preference != -9, it will extract the enriched cell-types within the cell-types identified with the naming preferences option.
  # as a final step, if an enriched cell-type is neuronal, it will complete a fisher's test of neuronal subtypes and replace the cell-marker with the most enrichriched neuronal subtype
  # Generally, this method seems to be biased to cell-types with a greater number of markers.
  # Args:
  # gene_lists is a named list of gene-symbols that is tested against -- when internal, it's extracted from topgenes
  # theSpecies: human or mouse
  # cell_marker_path: if local, path to Cell-Type marker rda files.
  # naming_preference: either -9 if there is no expected cell-type or one of the categories from get_naming_preference_options()
  # this is useful if you previously have an idea of which cell-type you were going to enrich
  # Returns:
  # MarkerSets: the gene set enrichment for each cell-type (to see what was the second/third most significant etc)
  # cellTypes: The top cell-type for each marker
  
  gmt_list <- "" # for CRAN dependencies
  cell_preference_final <- "" # for CRAN dependencies
  
  topGenes <- gene_lists
  marker_sets <- list()
  
  
  if(class(topGenes)== "character") {
    topGenes1 <- list()
    topGenes1[[1]] <- topGenes
    topGenes <- topGenes1
  }
  
  cellTypes <- c()
  
  for(i in 1:length(topGenes)) {
    #print(i)
    geneInput <- topGenes[[i]]
    if(theSpecies == "mouse") { # if it's a mouse
      #
      # Load in the cell-marker enrich stuff and switch names to fit cell-marker enrich
      #
      
      thefiles <- list.files(path = cell_marker_path, "cell_markers.rda")
      
      
      if(length(thefiles) == 0) {
        warning(paste0("Cell-marker databases are not present in ", cell_marker_path, " downloading and loading data."))
        #
        datafile <- "mouse_cell_markers.rda"
        metafile <- paste0(datafile)
        url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
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
        load(paste0(cell_marker_path,"/mouse_cell_markers.rda"))
        gmt_both <- gmt_list$gmt_both
        gmt_cellmarker <- gmt_list$gmt_cellmarker
        gmt_gobp <- gmt_list$gmt_gobp
        gmt_panglao <- gmt_list$gmt_panglao
        gmt_subtype <- gmt_list$gmt_subtype
      }

      outCellMarker <- cellmarker_enrich(topGenes[[i]], 0.05,gmt_cellmarker, fixed_length = 15000)
      # get the cell type from 'cellmarkers' dataset
      if(class(outCellMarker$name) == "factor") {
        outCellMarker$name <- tochr(outCellMarker$name)
      }
      outPanglao <- cellmarker_enrich(topGenes[[i]], 0.05, gmt_panglao, fixed_length = 15000)
      # get the cell type from 'panglao' dataset
      if(class(outPanglao$name) == "factor") {
        outPanglao$name <- tochr(outPanglao$name)
      }
      outGOMouse <- cellmarker_enrich(topGenes[[i]], 0.05, gmt_gobp, fixed_length = 22000)
      # get the cell type from gene ontology
      if(class(outGOMouse$name) == "factor") {
        outGOMouse$name <- tochr(outGOMouse$name)
      }
      
      marker_list <- list(CellMarker = outCellMarker, Panglao = outPanglao, Ontology = outGOMouse)
      # concatenate all of the CT markers from the appropriate dataset
      if(naming_preference != -9) {
        
        
        if(length(thefiles) == 0) {
          warning(paste0("Cell-marker preferences are not present in ", cell_marker_path, " downloading and loading data."))
          #
          datafile <- "cell_preferences_categorized.rda"
          metafile <- paste0(datafile)
          url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                        metafile, "?raw=true")
          destfile <- file.path(tempdir(), metafile)
          downloader::download(url, destfile = destfile, mode = "wb")
          load(destfile)
        } else {
          
          load(paste0(cell_marker_path,"/cell_preferences_categorized.rda"))
        }
        
        
        
        # if they have the prefered tissue or cell-type to pick from
        # prioritize those
        mypref <- cell_preference_final[[naming_preference]]
        outcell_preference <- which(outCellMarker$name %in% mypref)
        if(length(outcell_preference) > 0) {
          outCellMarker <- outCellMarker[outcell_preference,]
        }
        panglao_preference <- which(outPanglao$name %in% mypref)
        if(length(panglao_preference) > 0) {
          outPanglao <- outPanglao[panglao_preference,]
        }
        
      }
      marker_sets[[i]] <- marker_list 
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) == 0) {
        cellTypes[i] <- "unknown"
        # if none of that works out it's unknown
      }
      
      if(nrow(outGOMouse) == 0 & nrow(outGOMouse) == 0 & nrow(outGOMouse) > 0) {
        thenames <- sub( "\\%.*", "", outGOMouse$name)
        cellTypes[i] <- paste0("Ontology ", thenames[1])
        # define by top gene ontology if no other marker is present
        
      }
      
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) == 0) {
        thenames <- sub('.*\\:', '', outCellMarker$name)
        cellTypes[i] <- thenames[1]
        # define by cellmarker alone if none of the others are present
      }
      
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) > 0) {
        thenames <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- thenames[1]
        # define by panglao if it's the only dataset with a marker
      }
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) > 0) {
        thenamesCM <- sub('.*\\:', '', outCellMarker$name)
        thenamesP <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- paste0(thenamesCM[1], "_and_", thenamesP[1])
        # define by both if there was a signfiicant markerin both
      }
      
    }
    if(theSpecies == "human") {
      # this section is exactly the same as in mouse but generated for human
      # cell-marker lists
      thefiles <- list.files(path = cell_marker_path, "cell_markers.rda")
      
      
      if(length(thefiles) == 0) {
        warning(paste0("Cell-marker databases are not present in ", cell_marker_path, " downloading and loading data."))
        #
        datafile <- "human_cell_markers.rda"
        metafile <- paste0(datafile)
        url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
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
        load(paste0(cell_marker_path,"/human_cell_markers.rda"))
        gmt_both <- gmt_list$gmt_both
        gmt_cellmarker <- gmt_list$gmt_cellmarker
        gmt_gobp <- gmt_list$gmt_gobp
        gmt_panglao <- gmt_list$gmt_panglao
        gmt_subtype <- gmt_list$gmt_subtype
      }      
      outCellMarker <- cellmarker_enrich(topGenes[[i]], 0.05, gmt_cellmarker, fixed_length = 15000)
      if(class(outCellMarker$name) == "factor") {
        outCellMarker$name <- tochr(outCellMarker$name)
      }
      outPanglao <- cellmarker_enrich(topGenes[[i]], 0.05, gmt_panglao, fixed_length = 15000)
      if(class(outPanglao$name) == "factor") {
        outPanglao$name <- tochr(outPanglao$name)
      }
      outGoHuman <- cellmarker_enrich(topGenes[[i]], 0.05, gmt_gobp, fixed_length = 22000)
      if(class(outGoHuman$name) == "factor") {
        outGoHuman$name <- tochr(outGoHuman$name)
      }
      if(naming_preference != -9) {
        
        if(length(thefiles) == 0) {
          warning(paste0("Cell-marker preferences are not present in ", cell_marker_path, " downloading and loading data."))
          #
          datafile <- "cell_preferences_categorized.rda"
          metafile <- paste0(datafile)
          url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                        metafile, "?raw=true")
          destfile <- file.path(tempdir(), metafile)
          downloader::download(url, destfile = destfile, mode = "wb")
          load(destfile)
        } else {
          
          load(paste0(cell_marker_path,"/cell_preferences_categorized.rda"))
        }
        
        
        
        mypref <- cell_preference_final[[naming_preference]]
        outcell_preference <- which(outCellMarker$name %in% mypref)
        if(length(outcell_preference) > 0) {
          outCellMarker <- outCellMarker[outcell_preference,]
        }
        panglao_preference <- which(outPanglao$name %in% mypref)
        if(length(panglao_preference) > 0) {
          outPanglao <- outPanglao[panglao_preference,]
        }
        
      }
      marker_list <- list(CellMarker = outCellMarker, Panglao = outPanglao, Ontology = outGoHuman)
      marker_sets[[i]] <- marker_list 
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) == 0 & nrow(outGoHuman) == 0) {
        cellTypes[i] <- "unknown"
      }
      
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) == 0 & nrow(outGoHuman) > 0) {
        thenames <- sub( "\\%.*", "", outGoHuman$name)
        cellTypes[i] <- paste0("Ontology ", thenames[1])
        
      }
      
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) == 0) {
        thenames <- sub('.*\\:', '', outCellMarker$name)
        cellTypes[i] <- thenames[1]
      }
      
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) > 0) {
        thenames <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- thenames[1]
      }
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) > 0) {
        thenamesCM <- sub('.*\\:', '', outCellMarker$name)
        thenamesP <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- paste0(thenamesCM[1], "_and_", thenamesP[1])
      }
      
    }
  }

  names(marker_sets) <- cellTypes
  return(list(cellTypes = cellTypes, marker_sets = marker_sets))
}

