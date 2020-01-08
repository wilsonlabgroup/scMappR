#' Generate Heatmap
#'
#'  This function takes an inputted signature matrix as well as a list of genes and overlaps them.
#'  Then, if there is overlap, it prints a heatmap or barplot (depending on the number of overlapping genes).
#'  Then, for every cell-type, genes considered over-represented are saved in a list.
#'
#'
#' @rdname heatmap_generation
#' @name heatmap_generation
#'
#' @param genesIn A list of gene symbols (all caps) to have their cell type enrichment.
#' @param comp The name of the comparison.
#' @param cex The size of the genes in the column label for the heatmap.
#' @param rd_path The directory to RData files -- if they are not in this directory, then the files will be downloaded.
#' @param isMax If you are taking the single best CT marker (T/F) -- TRUE not reccomended.
#' @param isPval If the signature matrix is raw p-value (T/F) -- TRUE not reccomended.
#' @param cellTypes Colnames of the cell-types you will extract (passed to extract_genes_cell).
#' @param pVal The level of association a gene has within a cell type (passed to extract_genes_cell).
#' @param isBackground If the heatmap is from the entire signature matrix or just the inputted gene list (T/F). isBackground == TRUE is used for internal.
#' @param reference Path to signature matrix or the signature matrix itself.
#' @param which_species Species of gene symbols -- "human" or "mouse" .
#' @param toSave Allow scMappR to write files in the path directory (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#'
#' @return \code{heatmap_generation} A heatmap/barplot of p-value or odds-ratio of cell-type specific genes intersecting with the gene list. A list of genes that do/don't intersect with the signature matrix as well as a list of which cell-type these over-represented genes live in. \cr
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
#' # load in signature matrices
#' data(POA_example)
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' Signature <- POA_Rank_signature
#' rowname <- get_gene_symbol(Signature)
#' rownames(Signature) <- rowname$rowname
#' genes <- rownames(Signature)[1:100]
#' heatmap_test <- heatmap_generation(genesIn = genes, "scMappR_test",
#'                                    reference = Signature, which_species = "mouse")
#'
#' @export
#' 
heatmap_generation <- function(genesIn, comp,reference, cex = 0.8, rd_path = "~/scMappR/data", cellTypes = "ALL", pVal = 0.01, isPval=TRUE, isMax =FALSE,  isBackground = FALSE,  which_species = "human", toSave = FALSE, path = NULL) {
  # This function takes an inputted signature matrix as well as a list of genes and overlaps them. Then, if there is overlap, it prints a heatmap or barplot (depending on the number of overlapping genes)
  # Then, for every cell-type, genes considered over-represented are saved in a list
  
  # Args:
  # genesIn = a list of gene symbols (all caps) to have their cell type enrichment
  # comp = the name of the comparison
  # cex = the size of the genes in the column label for the heatmap... I know
  # The parameters below are passed into extract_genes_cell
  # cellTypes = which cell types you care about extracting (colnames)
  # pVal = the level of association a gene has within a cell type
  # isMax = true/false, just sort the gene into the cell type it's most strongly associated with
  # isBackground = if you're generating a map of an entire signature amtrix
  # refence = path to signature matrix if using an internal one, the actual signature matrix if custom
  # which species: "human" or "mouse" to designate the types of gene-symbols inputted 
  
  # Returns:
  # A heatmap/barplot of p-value or odds-ratio of cell-type specific genes intersecting with the gene list
  # a list of genes that do/don't intersect with the signature matrix as well as a list of which cell-type these over-represented genes live in.
  
  print(paste0("Assumes inputted gene symbols are of ", which_species, " origin."))
  if(class(genesIn) != "character") {
    stop("Inputted gene list is not character vector")
  }
  if(class(comp) != "character") {
    stop("comp must be of class caracter.")
  }
  if(class(cex) != "numeric") {
    stop("cex must be a numeric between 0 and 1")
  }
  if(cex < 0) {
    warning("cex < 0, setting to 0.001 (essentially invisible)") 
    cex <- 0.001
  }

  if(class(reference) != "data.frame" & class(reference) != "matrix") {
    stop("Reference must be of class data.frame or matrix.")
  }
  # cellTypes = "ALL", pVal = 0.01, isPval=TRUE, isMax =FALSE,  isBackground = FALSE,  which_species = "human", toSave = FALSE
  if(class(cellTypes) != "character") {
    stop("cellTypes should be of class character -- either column names of cell-types to include or 'ALL' -- ALL is reccomended.")
  }
  if(class(pVal) != "numeric") {
    stop("val must be of class numeric.")
  }
  if(!(any(is.logical(isMax), is.logical(isPval),is.logical(isBackground),is.logical(toSave)))) {
    stop("isMax, isBackground, toSave, and isPval must be of class logical (TRUE/FALSE).")
  }
  if(class(rd_path) != "character") {
    stop("rd_path must be of class character.")
  }
    
  if(!(which_species %in% c("human", "mouse"))) {
    stop("which_species must be either 'human' or 'mouse' (case sensitive).")
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  if(class(reference) == "data.frame") reference <- as.matrix(reference)
  if(class(reference) == "matrix") { # then you have a custom signature matrix or it's pre-loaded into R
    wilcoxon_rank_mat_t <- reference  
  }
    wilcoxon_rank_mat_t <- wilcoxon_rank_mat_t[!duplicated(rownames(wilcoxon_rank_mat_t)),]
    if(length(grep("-", rownames(wilcoxon_rank_mat_t))) / length(rownames(wilcoxon_rank_mat_t)) > 0.75) {  
      print("Detected signature matrix from scMappR catelogue", quote = FALSE)
      RN_2 <- get_gene_symbol(wilcoxon_rank_mat_t)
      
      rownames(wilcoxon_rank_mat_t) <- RN_2$rowname
      
    } else {
      RN_2 <- list(rowname = rownames(wilcoxon_rank_mat_t), species = which_species)
    }
    if(which_species != RN_2$species ) { # convert background species to your species (internal only)
      
      warning(paste0("Species in signature matrix, ",RN_2$species, " is not the same as the inputted gene species ", which_species,"."))
      warning(paste0( "Converting gene symbols in background from ",RN_2$species, " to ", which_species))
      
      
      
      thefiles <- list.files(path = rd_path, "bioMart_ortholog_human_mouse.rda")
      
      
      if(length(thefiles) == 0) {
        warning(paste0("Cell-marker databases are not present in ", rd_path, " downloading and loading data."))
        #
        datafile <- "bioMart_ortholog_human_mouse.rda"
        metafile <- paste0(datafile)
        url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                      metafile, "?raw=true")
        destfile <- file.path(tempdir(), metafile)
        downloader::download(url, destfile = destfile, mode = "wb")
        load(destfile)
        #
      } else {
        load(paste0(rd_path,"/bioMart_ortholog_human_mouse.rda"))
      
      
      rownames(bioMart_orthologs) <- bioMart_orthologs[,RN_2$species]
      RN <- rownames(wilcoxon_rank_mat_t)
      
      intersected_genes <- intersect(RN, rownames(bioMart_orthologs)) # gene sybmols in signature
      wilcoxon_gene1 <- wilcoxon_rank_mat_t[intersected_genes,]
      
      bioMart_orthologs1 <- bioMart_orthologs[intersected_genes,]
      rownames(wilcoxon_gene1) <- bioMart_orthologs1[,which_species] # replacing rownames with the species you want
      wilcoxon_rank_mat_t <- wilcoxon_gene1
      
    }
  }
  
  
  if(any(duplicated(colnames(wilcoxon_rank_mat_t)))) { # if cell-types are named the same (i.e. two clusters have the same tag then convert)
    print("Cell-types not uniquely named, appending tag to make all cell-types unique")
    colnames(wilcoxon_rank_mat_t) <- paste0(colnames(wilcoxon_rank_mat_t),"_tag",1:ncol(wilcoxon_rank_mat_t))
    
  }
  wilcoxon_rank_mat_t <- wilcoxon_rank_mat_t[!duplicated(rownames(wilcoxon_rank_mat_t)),] # just incase remove genes with the same name
  
  genesInter <- intersect(genesIn, rownames(wilcoxon_rank_mat_t)) # get genes with some DE
  whichGenesInter <- which(rownames(wilcoxon_rank_mat_t) %in% genesInter)
  genes_noInter <- genesIn[!(genesIn %in% genesInter)]
  geneHeat <- c()
  if(isBackground == TRUE) {
    print("Preferential expression in background set: ")
  }
  if(isBackground == FALSE) {
    print("Preferential expression in inputted gene list: ")
  }
  if(length(whichGenesInter) == 0) { # if no genes inputted are preferentially expressed
    print("No input genes are preferentially expressed")
    stop("No_downstram_analysis")
  }
  if(length(whichGenesInter) == 1) { # if only 1 gene is
    print("One input gene is preferentially expressed")
    if(toSave == TRUE) {
    grDevices::pdf(paste0(path,"/",comp,"_barplot.pdf"))
    graphics::barplot(wilcoxon_rank_mat_t[whichGenesInter,], las = 2, main = rownames(wilcoxon_rank_mat_t)[whichGenesInter])
    grDevices::dev.off()
    } else {
      warning("toSave = F and threfore plots are not allowed to be saved. I would reccomend allowing it to be true.")
    }
    stop("No_downstream_analysis")
  }
  if(length(whichGenesInter) > 1) { # if > 1 genes are
    print("At least one input gene is preferentially expressed")
    # make the heatmap
    if(toSave == TRUE) {
    myheatcol <- grDevices::colorRampPalette(c("lightblue", "white", "orange"))(256)
    grDevices::pdf(paste0(path,"/",comp,"_heatmap.pdf"))
    gplots::heatmap.2(wilcoxon_rank_mat_t[whichGenesInter,], Rowv = TRUE, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
    grDevices::dev.off()
    geneHeat <- wilcoxon_rank_mat_t[whichGenesInter,]
    preferences <- extract_genes_cell(geneHeat, cellTypes = cellTypes, val = pVal, isMax = isMax, isPvalue = isPval)

    save(preferences, file = paste0(path,"/",comp,"_preferences.RData"))
    } else {
      warning("toSave = F and scMappR is not allowed to print files or plots in your directories. For full functionality of the package, set to true.")
      geneHeat <- wilcoxon_rank_mat_t[whichGenesInter,]
      preferences <- extract_genes_cell(geneHeat, cellTypes = cellTypes, val = pVal, isMax = isMax, isPvalue = isPval)
      
    }
    
  } 
  
  
  return(list(genesIn = genesInter, genesNoIn = genes_noInter, geneHeat = geneHeat, preferences = preferences))
  # genesIn = genes with some preference, genesNoIn = genes with no preference 
  # geneHeat = matrix of preferences and p-values
  # preferences = genes mapped to their CT preference
}
