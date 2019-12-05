#' Count Matrix To Signature Matrix
#'
#' This function takes a list of count matrices, procces them, calls cell-types, and genreates signature matrices. much of the data is kept.
#' 
#' This function is a one line wrapper to process count matrices into a signature matrix.
#' It combines process from count, two methods of identifying cell-type identitt (gsva and fisher's test).
#' Then, it takes the output of cell-type markers and converts it into a signature matrix of P-value ranks and odds-ratios.
#' Along the way, it saves the Seurat object (if you choose), cell-type identites from gsva (it's own obect), and the signature matrices.
#' Cell-type marker outputs are also saved in the generes.RData list. Names of the generes objects and the signature matrices are kept.
#'
#' @rdname process_dgTMatrix_lists
#' @name process_dgTMatrix_lists
#'
#' @param dgTMatrix_list A list of matrices in the class of dgTMatrix object -- sparce object -- compatible with Seurat rownames should be of the same species for each.
#' @param name The name of the outputted signature matrices, cell-type preferences, and seurat objects if you choose to save them.
#' @param species_name Mouse or human symbols, -9 if internal as panglao objects have gene symbol and ensembl strapped together.
#' @param naming_preference For cell-type naming, see if cell-types given the inputted tissues are more likely to be named within one of the categories of get_naming_preference_options().
#' @param panglao_set If the inputted matrices are from panglao (i.e. if they're internal).
#' @param haveUMAP Save the umaps -- only use if the package is downloaded with pip.
#' @param saveSCObject Save the seurat object as an RData object (T/F).
#' @param internal Was this used as part of the internal processing of panglao datasets (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#' @param rda_path If saved, directory to where data from scMappR_data is downloaded.
#'
#' @return \code{process_from_count} A processed & integrated Seurat object that has been scaled and clustered. It can be returned as an internal object or also stored as an RData object if neccesary. \cr
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
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#'
#' @examples 
#' \donttest{
#' data(single_cell_process)
#' toProcess <- list(example = sm)
#' tst1 <- process_dgTMatrix_lists(toProcess, "testProcess", -9, 
#'                                 "eye",rda_path = "~/scMappR/data", TRUE)
#' }
#' @export
#' 

process_dgTMatrix_lists <- function(dgTMatrix_list, name, species_name, naming_preference,rda_path="",  panglao_set = FALSE ,haveUMAP = FALSE, saveSCObject = FALSE, internal = FALSE, toSave = FALSE) {
  
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
  sm <- dgTMatrix_list[[1]] # first thing we're going to do is figure out the species if it is currently unknown (mostly for internal)
  if(species_name == -9) {
    spec=get_gene_symbol(sm)
    species_name <- spec$species
  } 
  
  pbmc <- process_from_count(countmat_list = dgTMatrix_list, name = name, theSpecies = species_name, panglao_set = panglao_set, haveUmap = haveUMAP, saveALL  = saveSCObject, toSave=toSave)
  # process from the count matrices to the Seurat object -- see process_from_count for details
  print(class(pbmc))
  #print(head(pbmc))
  gsva_cellIdentity_out <- gsva_cellIdentify(pbmc, theSpecies = species_name, naming_preference = naming_preference, rda_path = rda_path, toSave=toSave)
  # identify cell-type identify using the gsva method
  if(toSave == TRUE) {
    save(gsva_cellIdentity_out, file = paste0(name, "gsva_cellname.Rdata"))
      
  } else {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.")
    print("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.", quote = FALSE)
    
  }
  generes <- seurat_to_generes(pbmc = pbmc)
  # identify cell-type marker for every cluster identified
  gene_out <- generes_to_heatmap(generes, species = species_name, naming_preference = naming_preference, internal = internal, rda_path = rda_path)
  # generate signature matrices and identify cell-types with the fisher's test method
  wilcoxon_rank_mat_t <- gene_out$pVal
  wilcoxon_rank_mat_or <- gene_out$OR
  names(generes) <- colnames(wilcoxon_rank_mat_t) <- colnames(wilcoxon_rank_mat_or) <- gene_out$cellname
  # save cell-type markers and signature matrices
  if(toSave == TRUE) {
    save(generes, file = paste0(name, "_generes.Rdata"))
    save(wilcoxon_rank_mat_t, file= paste0(name, "_pval_heatmap.Rdata"))
    save(wilcoxon_rank_mat_or, file= paste0(name, "_or_heatmap.Rdata"))
    l <- list(wilcoxon_rank_mat_t = wilcoxon_rank_mat_t, wilcoxon_rank_mat_or = wilcoxon_rank_mat_or,generes=generes)
  } else {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.")
    print("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.", quote = FALSE)
    
    l <- list(wilcoxon_rank_mat_t = wilcoxon_rank_mat_t, wilcoxon_rank_mat_or = wilcoxon_rank_mat_or,generes=generes)
    
  }
  
  return(l)
}

