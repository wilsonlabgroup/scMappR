#' Count Matrix To Seurat Object
#'
#' This function processes a list of count matrices (same species/gene symbols in each list) and converts them to a Seurat object
#' 
#' This function takes a list of count matrices and returns a seurat object of the count matrices integrated using Seurat V3 and the interation anchors
#' Different options are used for if the function is internal for PanglaoDB dataset reprocessing or being used for a custom set of count matrices.
#' For larger scRNA-seq datasets (~20k + cells), it is likely that this function will be required to run on an hpc.
#'
#' @rdname process_from_count
#' @name process_from_count
#'
#' @param countmat_list A list of count matrices that will be be integrated using the integration-anchors features they should have the same rownames.
#' @param name The output of the normalzied and fused Suerat object if you choose to keep it.
#' @param theSpecies Gene symbols for human, mouse, or -9 if internal. If your species is not human or mouse gene symbols, make sure that you have "MT-" before your mitochondrial gene names then pick "human".
#' @param haveUmap Write a UMAP (T/F).
#' @param saveALL Save the Seurat object generated (T/F).
#' @param panglao_set If the function is being used from internal (T/F).
#' @param toSave Allows scMappR to print files and make directories locally (T/F).
#' @param use_sctransform If you should use scRNAsform or the normalize/variablefeatures/scaledata pipeline (T/F).
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
#' @importFrom gprofiler2 gost
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#'
#' @examples 
#' 
#' data(sm)
#' toProcess <- list(example = sm)
#' tst1 <- process_from_count(toProcess, "testProcess")
#' 
#' @export
#' 
process_from_count <- function(countmat_list, name, theSpecies = -9, haveUmap = FALSE, saveALL = FALSE, panglao_set = FALSE, toSave = FALSE, use_sctransform = FALSE) {
  # This function takes a list of count matrices and returns a seurat object of the count matrices integrated using Seurat V3 and the interation anchors
  # Different options are used for if the function is internal for PanglaoDB dataset reprocessing or being used for a custom set of count matrices.
  # For larger scRNA-seq datasets (~20k + cells), it is likely that this function will be required to run on an hpc.
  
  # Args:
  # countmat_list: a list of count matrices that will be be integrated using the integration-anchors features
  # they should have the same rownames
  # name = the output of the normalzied and fused Suerat object if you choose to keep it.
  # theSpecies = gene symbols for human, mouse, or -9 if internal. If your species is not human or mouse gene symbols
  #  then insure that your mitochondiral genes contain "MT-" at the beginning on the gene and pick human.
  # haveUmap =  write a Umap T/F
  # saveALL = save the Seurat object generated T/F
  # Panglao_set = If the function is being used from internal T/F
  # use_sctransform = If you should use scRNAsform or the normalize/variablefeatures/scaledata pipeline (T/F)
  # Returns:
  # A processed & integrated Seurat object that has been scaled and clustered. It can be returned as an internal object or 
  # also stored as an RData object if neccesary.
  
  if(!is.character(name)) {
    stop("Name is not a character for your outputs, please change the parameter and try again.")
  }
    
  
  if(class(countmat_list) == "dgCMatrix" | class(countmat_list) == "matrix") {
    print("'dgTMatrix_list' is of class dgCMatrix or matrix, converting to a named list.", quote = F)
    countmat_list <- list(name = countmat_list)
    names(countmat_list) <- name
  }
  if(is.null(names(countmat_list))) {
    warning("List has no names, adding names")
    names(countmat_list) <- paste0(name,"_",1:length(countmat_list))
  }
  
  SRA_in <- countmat_list
  
  shrt <- names(SRA_in)
  name <- name
  each_sra <- list()
  count <- 1
  for(f in 1:length(SRA_in)) { # for each count matrix
    
    sm <- SRA_in[[f]]
    sm <- sm[!duplicated(rownames(sm)),]
    ####################################
    ####################################
    
    if(theSpecies == -9 | panglao_set == TRUE) { # if internal and we want to shave ENSEMBL symbols from rownames
      sym <- get_gene_symbol(sm)
      RN_2 <- sym$rowname
      theSpecies <- sym$species
      rownames(sm) <- RN_2
      
    }
    RN_2 <- rownames(sm)
    if(theSpecies == "mouse") { # test to see if mitochondrial genes are designated
      num_MT <- grep("mt-", RN_2)
    }
    if(theSpecies == "human") {
      num_MT <- grep("MT-",RN_2)
    }
    
    # the code below will check to see if there are mt genes in the correct format.
    # If they are in the right format then continue on, 
    # otherwise adjust the gene names so that mito genes are detected
    mt.genes_m <- c("Tf", "Rnr1","Tv","Rnr2","Tl1","Nd1","Ti","Tq","Tm","Nd2","Tw","Ta","Tn","Tc","Ty","Co1","Ts1","Td","Co2","Tk","Atp8","Atp6","Co3","Tg","Nd3","Tr","Nd4l","Nd4","Ts2","Tl2","Nd5","Nd6","Te","Cytb","Tt","Tp")
    print(length(num_MT))
    if(length(num_MT) == 0 & theSpecies =="human") {
      # convert gene names if there are none with the mitochondiral designation.
      mt.genes <- toupper(mt.genes_m)
      mito.genes <- which(RN_2 %in% mt.genes)
      RN_2[mito.genes] <- paste0("MT-", RN_2[mito.genes])
      
      mito.genes <- RN_2[mito.genes]
    }
    
    if(length(num_MT) == 0 & theSpecies =="mouse") {
      # convert gene names mouse if there are none with the mitochondrial designation
      mt.genes <- mt.genes_m
      mito.genes <- which(RN_2 %in% mt.genes)
      mito.genes <- RN_2[mito.genes]
      RN_2[mito.genes] <- paste0("mt-", RN_2[mito.genes])
      
    }
    rownames(sm) <- RN_2
    
    pbmc <- Seurat::CreateSeuratObject(sm, min.cells = 3, min.features = 0, project = shrt[count])
    if(theSpecies == "human") {
      pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt.adj")
    }
    if(theSpecies == "mouse") {
      pbmc <- Seurat::PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "percent.mt.adj")
    }
    
    ## Remove cells with >2 standard deviations of MT contamination given the dataset.
    mean <- mean(pbmc$percent.mt.adj)
    x<- stats::sd(pbmc$percent.mt.adj)
    toremove <- mean + (2*x)
    
    if(use_sctransform == FALSE) { # alternative to using scTransform, use the traditional NormalizeData, FindVariableFeatures, and SaleData parameters. 
                                   # This should be faster and cost less memory
      
    if(toremove > 0) {
      pbmc <- pbmc[,which(pbmc$percent.mt.adj < toremove)]
      pbmc <- Seurat::NormalizeData(object = pbmc, normalization.method = "LogNormalize",
                                    scale.factor = 10000)
      pbmc <- Seurat::FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
      
      pbmc <- Seurat::ScaleData(pbmc, vars.to.regress = "percent.mt.adj")
    } else {
      
      pbmc <- Seurat::NormalizeData(object = pbmc, normalization.method = "LogNormalize",
                                    scale.factor = 10000)
      pbmc <- Seurat::FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                                           x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)
      
      pbmc <- Seurat::ScaleData(pbmc)
      
    }

    
    
    } else {
      
      if(toremove > 0){
        pbmc <- pbmc[,which(pbmc$percent.mt.adj < toremove)]
        pbmc <- Seurat::SCTransform(pbmc, vars.to.regress = "percent.mt.adj", verbose = FALSE)
      } else {
        warning("Dataset has no MT contamination -- this is highly unlikely. Please make sure that MT genes are designated with 'MT-' for human and 'mt-' for mouse. Normalizing on depth and not % mitochondria.")
        pbmc <- Seurat::SCTransform(pbmc, verbose = FALSE)
        
      }
      
    }
    #########
    #
    # Cannot process with scTransform because it doesn't allow for integration yet -- may be able to update
    
    
    
    # process with sctransformed (most updated)
    
    
    
    each_sra[[count]] <- pbmc
    count <- count + 1
    print(count)
    
  }
  names(each_sra) <- shrt
  
  if(length(SRA_in) > 1) {
    # If there is more than one count matrix in the list, then integrate it using the 
    # integration anchors feature using default parameters
    
    object.features <- Seurat::SelectIntegrationFeatures(object.list = each_sra, nfeatures = 3000)
    object.list <- Seurat::PrepSCTIntegration(object.list = each_sra, anchor.features = object.features, 
                                              verbose = FALSE)
    immune.anchors <- Seurat::FindIntegrationAnchors(object.list = object.list, anchor.features = object.features, normalization.method = "SCT")
    immune.combined <- Seurat::IntegrateData(anchorset = immune.anchors, dims = 1:20, normalization.method = "SCT")
    
    Seurat::DefaultAssay(immune.combined) <- "integrated"
    pbmc <- immune.combined
  }
  
  #PCA, UMAP, and Clustering of integrated dataset
  pbmc <- Seurat::RunPCA(object = pbmc, verbose = FALSE)
  if(haveUmap == TRUE) {
    pbmc <- Seurat::RunUMAP(object = pbmc, dims = 1:20, verbose = FALSE)
  }
  pbmc <- Seurat::FindNeighbors(object = pbmc, dims = 1:20, verbose = FALSE)
  pbmc <- Seurat::FindClusters(object = pbmc, verbose = FALSE)
  if(saveALL == TRUE & toSave == TRUE) {
    # Save the seurat object before scaling
    save(pbmc, file = paste0(name, "_custom.Rdata"))
  }
  pbmc <- try(Seurat::ScaleData(object = pbmc)) # scale data
  if(class(pbmc) == "try-error") {
    stop("Data scaling did not finish, this can often be due to a memory error (as of November 2019). 
        the Seurat object up to this point has been saved.
        SCTransform can lead to this memory issue. consider the 'process from count no sctransform'
         from the scMappR github if you can't get it to work by chaging memory options" )
    
  }
  if(saveALL == TRUE & toSave == TRUE) {
    
    save(pbmc, file = paste0(name, "_custom.Rdata"))
  } else {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended.")
  }
  return(pbmc)
}

