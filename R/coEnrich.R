#' Identify co-expressed cell-types
#'
#' This function identifies genes with similar cell-type markers and if those markers are driving enrichment.
#'
#' This function takes significantly enriched cell-types from the single CT_enrich before testing to see if the genes driving their enrichment are overlapping to a significant proportion using Fisher's exact test.
#' To save computational time and to not complete this with an incredible number of permutations, scMappR stops at overlapping 5 cell-types.
#'
#'
#' @rdname coEnrich
#' @name coEnrich
#'
#' @param sig A The number of combinations of significant cell-types to enrich.
#' @param gene_list_heatmap Signature matrix of inputted genes in heatmap and the cell-type preferences -- output of heatmap generation.
#' @param background_heatmap Signature matrix of background matrix in heatmap and cell-type preferences -- output of heatmap generation.
#' @param study_name Name of the outputted table.
#' @param outDir Name of the directory this table will be printed in.
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#'
#' @return \code{coEnrich} Enrichment of cell-types that are expressed by the same genes, up to 4 sets of cell-types. \cr
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
#'
#' @examples 
#' \donttest{
#' # load in signature matrices
#' data(POA_example)
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' sig <- get_gene_symbol(POA_Rank_signature)
#' Signature <- POA_Rank_signature
#' rownames(Signature) <- sig$rowname
#' genes <- rownames(Signature)[1:60]
#' heatmap_test <- tissue_scMappR_custom( genes, signature_matrix = Signature,
#' output_directory =  "scMappR_test", toSave = FALSE)
#' group_preferences <- heatmap_test$group_celltype_preferences
#' }
#' @export
#' 
coEnrich <- function(sig, gene_list_heatmap, background_heatmap, study_name, outDir, toSave = FALSE, path = NULL) {
  # Internal
  # this function takes significantly enriched cell-types from the single CT enrich before testing to see
  # if the genes driving their enrichment are overlapping
  # to save computational time and to not complete this with an incredible number of permutations, scMappR stops at overlapping 5 cell-types
  # Args:
  # sig = preferences of significant functions
  # gene_list_heatmap = signature matrix (p-value) of your subsetted gene list
  # background_heatmap = signature matrix (p-value) of entire tissue
  # study_name = the name of the co-enrich table to be outputted
  # output_directory = the name of the directory where co-enrich will be completed
  # Returns:
  # Enrichment of cell-types that are expressed by the same genes, from 2-5 sets of cell-types.
  
  if((any(is.null(sig), is.null(gene_list_heatmap), is.null(background_heatmap), is.null(study_name), is.null(outDir)))[1]) {
    stop("One of the arguments is NULL suggesting that this function is not being run internally.")
  } 
  if(!is.logical(toSave)) {
    stop("toSave is not a logical object (TRUE/FALSE)")
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  if((nrow(sig) > 5)[1]) {
    sig <- sig[1:5,]
  }
  l <- nrow(sig)

  multi_comps <- c()
  for(y in 2:l) {
    # for combinations of 2-# enriched cell types (max= 5)
    if(y <= l) {
      comps <- utils::combn(sig$cell_type, y)
      # find combinations of cell-types
      co_up <- function(x) return(length(x[x>=1])==y)
      #score genes significantly enriched in both cell types = 1
      message(y)
      for(j in 1:ncol(comps)) {
        # for that comparison
        thecomps <- toupper(comps[,j]) # take comparison
        geneList_comb <- gene_list_heatmap$geneHeat # get gene list signature matrix
        colnames(geneList_comb) <- toupper(colnames(geneList_comb))
        geneList_comb1 <- geneList_comb[,which(colnames(geneList_comb) %in% thecomps)] # extract cell-types to be cenriched
        background_comb <- background_heatmap$geneHeat # take tissue signature matrix
        colnames(background_comb) <- toupper(colnames(background_comb))
        
        inter <- intersect(rownames(background_heatmap), rownames(geneList_comb))
        
        background_comb1 <- background_comb[!(rownames(background_comb) %in% inter), which(colnames(background_comb) %in% thecomps)] # remove enriched CT
  
        bin_aa <- apply(geneList_comb1,1,co_up) # see who is co-enriched in gene list
        name_in <- names(bin_aa)[bin_aa == T]
        aa <- sum(bin_aa) # number of co-enriched genes in list
        ab <- nrow(geneList_comb1) - aa # number non-co-enriched genes not in list
        ba <- sum(apply(background_comb1,1,co_up)) # number of co-enriched not in list
        bb <- nrow(background_comb1) - ba # number of non-coenriched genes not in list
        m <- matrix(c(aa, ba, ab, bb), nrow = 2) # fishers exact test
        fisherTest <- stats::fisher.test(m)
        OR <- fisherTest$estimate
        p <- fisherTest$p.value
        name <- paste0(thecomps,collapse=":" )
        row <- c(name,p,OR, paste0(name_in, collapse = ","))
        multi_comps <- rbind(multi_comps, row)
      }
    }
  }
  colnames(multi_comps) <- c("cell_types", "p_val", "OR", "genes")
  multi_comps <- as.data.frame(multi_comps)
  multi_comps$p_val <- toNum(multi_comps$p_val)
  multi_comps$pFDR <- p.adjust(multi_comps$p_val, "fdr")
  if(toSave == TRUE) {
  utils::write.table(multi_comps, file = paste0(path,"/",outDir, "/",study_name, "cell_co_preferences.tsv"), quote = FALSE, row.names = FALSE, col.names = T, sep = "\t")
  } else {
    warning("You are not allowing scMappR to save files. We strongly reccomend you switch toSave = TRUE")
  }
  
  return(multi_comps)
}

