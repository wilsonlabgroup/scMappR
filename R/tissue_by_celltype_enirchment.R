#' tissue_by_celltype_enrichment
#'
#' This function uses a Fisher's-exact-test to rank gene-set enrichment.
#' 
#' Complete a Fisher's exact test of an input list of genes against one of the two curated tissue by cell-type marker datasets from scMappR.
#'
#' @rdname tissue_by_celltype_enrichment
#' @name tissue_by_celltype_enrichment
#'
#' @param gene_list A character vector of gene symbols with the same designation (e.g. mouse symbol - mouse, human symbol - human) as the gene set database.
#' @param species Species of cell-type marker to use ('human' or 'mouse').
#' @param p_thresh The Fisher's test cut-off for a cell-marker to be enriched.
#' @param rda_path Path to a .rda file containing an object called "gmt". Either human or mouse cell-type markers split by experiment. If the correct file isn't present they will be downloaded from https://github.com/wilsonlabgroup/scMappR_Data.
#' @param isect_size Number of genes in your list and the cell-type.
#' @param return_gmt Return .gmt file -- reccomended if downloading from online as it may have updated (T/F).
#' @param name Name of the pdf to be printed.
#' 
#' @return List with the following elements:
#' \item{enriched}{Data frame of enriched cell-types from tissues.}
#' \item{gmt}{Cell-markers in enriched cell-types from tissues.}
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text geom_barplot theme_classic xlab ylab scale_fill_manual element_line
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
#' data(POA_example)
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' Signature <- POA_Rank_signature
#' rowname <- get_gene_symbol(Signature)
#' rownames(Signature) <- rowname$rowname
#' genes <- rownames(Signature)[1:100]
#' 
#' enriched <- tissue_by_celltype_enrichment(gene_list = genes, 
#' species = "mouse",p_thresh = 0.05, isect_size = 3)
#' 
#' }
#'  
#' @export
#' 

tissue_by_celltype_enrichment <- function(gene_list, species, name = "CT_Tissue_example", p_thresh = 0.05, rda_path = "",  isect_size = 3, return_gmt= FALSE) {
  gmt <- "" # no visible binding 
  if(is.null(species)) stop("please select 'human' or 'mouse' as a species.")
  
  if(!is.character(gene_list)) {
    stop("Gene list is not a character vector.")
  }
  if(length(gene_list) == 0) {
    stop("length of the gene list is 0, please check and try again")
  }
  if(length(gene_list) < 5) {
    warning("Fewer than 5 genes, cell-type enrichment may be challenging to identify.")
  }
  species_in <- species %in% c("human", "mouse")
  if(species_in[1] == FALSE) { # making sure species is correct
    stop("species does not equal human or mouse (case sensitive)")
  }
  
  if(!is.character(name)) {
    stop("name must be of class character.")
  }
  
  if(!is.character(rda_path)) {
    stop("rda_path must be of class character.")
  }
  
  if(!is.numeric(p_thresh) ) {
    stop("p_thresh must be of class numeric.")
  }
  
  if(!is.numeric(isect_size) ) {
    stop("isect_size must be of class numeric.")
  }
  if(all(is.logical(return_gmt))[1] == FALSE) {
    stop("return_gmt must be of class logical.")
  }
  

  if(species == "human") { # downloading human CT markers
  #
  thefiles <- list.files(path = rda_path, "human_tissue_celltype_scMappR.rda")
  if(length(thefiles) == 0) {  
    warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
    warning("Consider setting return_gmt to true as the gmt file may have updated")
    datafile <- "human_tissue_celltype_scMappR.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
  } else { # load cell-type markers locally
    load(paste0(rda_path,"/human_tissue_celltype_scMappR.rda"))
    
  }
  #
  }
  if(species == "mouse") { # download mouse CT markers
    #
    thefiles <- list.files(path = rda_path, "mouse_tissue_celltype_scMappR.rda")
    if(length(thefiles) == 0) {
    warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
    warning("Consider setting return_gmt to true as the gmt file may have updated")
      
    datafile <- "mouse_tissue_celltype_scMappR.rda"
    metafile <- paste0(datafile)
    url <- paste0("https://github.com/wilsonlabgroup/scMappR_Data/blob/master/", 
                  metafile, "?raw=true")
    destfile <- file.path(tempdir(), metafile)
    downloader::download(url, destfile = destfile, mode = "wb")
    load(destfile)
    } else { # you already had it, load locally
      load(paste0(rda_path,"/mouse_tissue_celltype_scMappR.rda"))
      
    }
  }
  
  enriched <- cellmarker_enrich(gene_list = gene_list, p_thresh = p_thresh, gmt = gmt, isect_size = isect_size)
  if(nrow(enriched) == 0) {
    warning("No cell-type markers were enriched, consider adjusting p_thresh and isect_size as well as making sure that the gene symbols in your input gene list matches the gmt.")
  }
  enriched$term_name <- tochr(enriched$name)
  enriched$p_value <- enriched$fdr

  if(return_gmt == TRUE) {
    l <- list(enriched = enriched, gmt = gmt) # returning the enrichment of the gene list and background gmt as it updates 
    return(l)
  } else {
    return(enriched)
  }
}

