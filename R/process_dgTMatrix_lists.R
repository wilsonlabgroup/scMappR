#' Count Matrix To Signature Matrix
#'
#' This function takes a list of count matrices, processes them, calls cell-types, and generates signature matrices. 
#' 
#' This function is a one line wrapper to process count matrices into a signature matrix.
#' It combines process_from_count, two methods of identifying cell-type identities (GSVA and Fisher's test).
#' Then, it takes the output of cell-type markers and converts it into a signature matrix of p-value ranks and odds ratios.
#' Along the way, it saves the Seurat object (if chosen with saveSCObject), cell-type identites from GSVA (its own obect), and the signature matrices.
#' Cell-type marker outputs are also saved in the generes.RData list. Names of the generes objects and the signature matrices are kept.
#'
#' @rdname process_dgTMatrix_lists
#' @name process_dgTMatrix_lists
#'
#' @param dgTMatrix_list A list of matrices in the class of dgTMatrix object -- sparce object -- compatible with Seurat rownames should be of the same species for each.
#' @param name The name of the outputted signature matrices, cell-type preferences, and seurat objects if you choose to save them.
#' @param species_name Mouse or human symbols, -9 if internal as panglao objects have gene symbol and ensembl strapped together.
#' @param naming_preference For cell-type naming, see if cell-types given the inputted tissues are more likely to be named within one of the categories of get_naming_preference_options().
#' @param panglao_set If the inputted matrices are from Panglao (i.e. if they're internal).
#' @param haveUMAP Save the UMAPs -- only use if the package is downloaded with pip.
#' @param saveSCObject Save the Seurat object as an RData object (T/F).
#' @param use_sctransform If you should use sctransform or the Normalize/VariableFeatures/ScaleData pipeline (T/F).
#' @param test_ctname statistical test for calling CT markers -- must be in Seurat
#' @param internal Was this used as part of the internal processing of Panglao datasets (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#' @param rda_path If saved, directory to where data from scMappR_data is downloaded.
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#' @param genes_integrate The number of genes to include in the integration anchors feature when combining datasets -- passed into process_from_count
#' @param genes_include TRUE or FALSE -- include 2000 genes in signature matrix or all matrix.
#' 
#' @return List with the following elements:
#' \item{wilcoxon_rank_mat_t}{A dataframe containing the signature matrix of ranks (-log10(Padj) * sign(fold-change)).}
#' \item{wilcoxon_rank_mat_or}{A dataframe containing the signature matrix of odds-ratios.}
#' \item{generes}{All cell-type markers for each cell-type with p-value and fold changes.} 
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
#' data(sm)
#' toProcess <- list(example = sm)
#' tst1 <- process_dgTMatrix_lists(toProcess, name = "testPropcess", species_name = -9,
#'  naming_preference = "eye", rda_path = "")
#' 
#' }
#' @export
#' 

process_dgTMatrix_lists <- function(dgTMatrix_list, name, species_name, naming_preference = -9,rda_path="",  panglao_set = FALSE ,haveUMAP = FALSE, saveSCObject = FALSE, internal = FALSE, toSave = FALSE, path = NULL, use_sctransform = FALSE, test_ctname = "wilcox", genes_integrate = 2000, genes_include = FALSE) {
  
  # This function is a one line wrapper to process count matrices into a signature matrix
  # It combines process from count, two methods of identifying cell-type identitt (gsva and fisher's test)
  # then it takes the output of cell-type markers and converts it into a signature matrix of P-value ranks and odds-ratios
  # Along the way it saves the Seurat object (if you choose), cell-type identites from gsva (it's own obect), and the signature matrices
  # cell-type marker outputs are also saved in the generes.RData list. names of the generes objects and the signature matrices are kept.
  # Args:
  # dgTMatrix_list: a list of matrices in the class of dgTMatrix object -- sparce object -- compatible with Seurat rownames should be of the same species for each
  # name: the name of the outputted signature matrices, cell-type preferences, and seurat objects if you choose to save them
  # species: mouse or human symbols, -9 if internal as panglao objects have gene symbol and ensembl strapped together
  # naming_preference: for cell-type naming, see if cell-types given the inputted tissues are more likely to be named within one of the categories of get_naming_preference_options()
  # panglao_set: If the inputted matrices are from panglao (i.e. if they're internal)
  # have_umap: save the umaps -- only use if the package is downloaded with pip
  # saveSCobject: save the seurat object as an RData object (T/F)
  # internal: was this used as part of the internal processing of panglao datasets (T/F)
  # Returns:
  # Seurat object in Rdata object (optional)
  # cell-type markers from gsva in RData object
  # list of cell-type markers and summary statistics (P-val-adust and log2 fold-change)
  # these are also named by the fisher's exact test method
  # Signature matrices using a -log10(Padj) and fold-changes # in Rdata file
  # It returns the signature matrix of P-values as an R object. 
  
  if(!is.character(name)) {
    stop("Name is not a character for your outputs, please change the parameter and try again.")
  }
  
  dgTMatrix_list_class1 <- class(dgTMatrix_list)[1] %in% c("dgCMatrix", "matrix", "list")
  if(dgTMatrix_list_class1[1] == FALSE) {
    stop("'dgTMatrix_list' is not dgCMatrix, matrix, or list. Please input data in appropriate class.")
  }
  
  dgTMatrix_list_class2 <- class(dgTMatrix_list)[1] %in% c("dgCMatrix", "matrix")
  
  if(dgTMatrix_list_class2[1]) {
    message("'dgTMatrix_list' is of class dgCMatrix or matrix, converting to a named list.", quote = F)
    dgTMatrix_list <- list(name = dgTMatrix_list)
    names(dgTMatrix_list) <- name
  }
  if(is.null(names(dgTMatrix_list))) {
    warning("List has no names, adding names")
    names(dgTMatrix_list) <- paste0(name,"_",1:length(dgTMatrix_list))
  }
  
  sm <- dgTMatrix_list[[1]] # first thing we're going to do is figure out the species if it is currently unknown (mostly for internal)
  
  if(!(species_name %in% c("human", "mouse"))) {
    if(species_name != -9) {
      stop("species_name is not 'human' 'mouse' or '-9' (case sensitive), please try again with this filled.")
    }
  }
  
  if(!is.character(rda_path)) {
    stop("rda_path must be of class character.")
  }
  
  # panglao_set = FALSE ,haveUMAP = FALSE, saveSCObject = FALSE, internal = FALSE, toSave = FALSE, use_sctransform = FALSE)

  if(all(is.logical(panglao_set),is.logical(haveUMAP),is.logical(saveSCObject),is.logical(internal),is.logical(panglao_set),is.logical(toSave),is.logical(use_sctransform)) == FALSE) {
    stop("panglao_set, haveUMAP, saveSCObject, internal, toSave, and use_sctransform must all be of class logical.")
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  
    
  if(species_name == -9) {
    spec=get_gene_symbol(sm)
    species_name <- spec$species
    for(z in 1:length(dgTMatrix_list)) {
      sm_z <- dgTMatrix_list[[z]]
      Row_name <- get_gene_symbol(sm_z)
      row_name_use <- Row_name$rowname
      row_name_use_1 <- row_name_use[!duplicated(row_name_use)]
      sm_z <- sm_z[!duplicated(row_name_use),]
      rownames(sm_z) <- row_name_use_1
      dgTMatrix_list[[z]] <- sm_z
    }
  } 
  naming_preferences <- c("brain", "epithelial", "endothelial", "blood", "connective","eye", "epidermis", "Digestive", "Immune", "pancreas", "liver", "reproductive", "kidney", "respiratory") 
  if(!naming_preference %in% naming_preferences) {
    if(naming_preference != -9) {
      message("Naming preference options")
      message(naming_preferences)
      stop("Naming preferences not in options (case sensitive) and isn't a non-choice (-9), please try again.")
    }
  }
  
  pbmc <- process_from_count(countmat_list = dgTMatrix_list, name = name, theSpecies = species_name, panglao_set = panglao_set, haveUmap = haveUMAP, saveALL  = saveSCObject, toSave=toSave, use_sctransform = use_sctransform, path = path, genes_integrate = genes_integrate, genes_include= genes_include)
  # process from the count matrices to the Seurat object -- see process_from_count for details
  message(class(pbmc))
  #message(head(pbmc))
  message(naming_preference)
  message(class(naming_preference))
  gsva_cellIdentity_out <- gsva_cellIdentify(pbmc, theSpecies = species_name, naming_preference = naming_preference, rda_path = rda_path, toSave=toSave)
  # identify cell-type identify using the gsva method
  if(toSave == TRUE) {
    save(gsva_cellIdentity_out, file = paste0(path, "/", name, "_gsva_cellname_and_avg_expression.Rdata"))
      
  } else {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.")
    message("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.", quote = FALSE)
    
  }
  generes <- seurat_to_generes(pbmc = pbmc, test = test_ctname)
  # identify cell-type marker for every cluster identified
  gene_out <- generes_to_heatmap(generes, species = species_name, naming_preference = naming_preference, internal = internal, rda_path = rda_path)
  # generate signature matrices and identify cell-types with the fisher's test method
  wilcoxon_rank_mat_t <- gene_out$pVal
  wilcoxon_rank_mat_or <- gene_out$OR
  names(generes) <- colnames(wilcoxon_rank_mat_t) <- colnames(wilcoxon_rank_mat_or) <- gene_out$cellname
  # save cell-type markers and signature matrices
  
  # Combine all cell-type makers into a table
  
    # Primary cell-type names
    OR_names <- colnames(wilcoxon_rank_mat_or)
    OR_split <- strsplit(OR_names,"_and_")
    for(p in 1:length(OR_split)) {
      if(length(OR_split[[p]]) == 1) {
        OR_split[[p]][2] <- OR_split[[p]][1]
      }
    }
    OR_split_mat <- do.call("rbind",OR_split)
    
    # GSVA identity
    cellmarker <- unlist(unname(gsva_cellIdentity_out$cellMarker))
    panglao <- unlist(unname(gsva_cellIdentity_out$panglao))
  
    # topGenes
    topGenes <- scMappR::topgenes_extract(generes)
    
    colPaste <- function(x) return(paste0(x,sep=",",collapse=""))
    topGenes_paste <- unname(unlist(lapply(topGenes, colPaste) ))
    cluster_num <- 1:length(topGenes_paste) - 1
    name_together <-  cbind(cluster_num, OR_split_mat,cellmarker,panglao,topGenes_paste)
    colnames(name_together) <- c("clusterNum", "CellMarkerFisher", "PanglaoFisher", "CellMarkerGSVA", "PanglaoGSVA", "TopGenes")
    
  ##
  
  if(toSave == TRUE) {
    save(generes, file = paste0(path,"/", name, "_generes.Rdata"))
    save(wilcoxon_rank_mat_t, file= paste0(path,"/", name, "_pval_heatmap.Rdata"))
    save(wilcoxon_rank_mat_or, file= paste0(path,"/", name, "_or_heatmap.Rdata"))
    save(name_together, file = paste0(path,"/",name,"_cellLabelInfo.Rdata"))
    l <- list(wilcoxon_rank_mat_t = wilcoxon_rank_mat_t, wilcoxon_rank_mat_or = wilcoxon_rank_mat_or,generes=generes, cellLabel = name_together)
  } else {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.")
    message("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.", quote = FALSE)
    
    l <- list(wilcoxon_rank_mat_t = wilcoxon_rank_mat_t, wilcoxon_rank_mat_or = wilcoxon_rank_mat_or,generes=generes, cellLabel = name_together)
    
  }
  
  return(l)
}

