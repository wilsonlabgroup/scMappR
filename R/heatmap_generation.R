#' Generate Heatmap
#'
#'  This function takes an inputted signature matrix as well as a list of genes and overlaps them.
#'  Then, if there is overlap, it prints a heatmap or barplot (depending on the number of overlapping genes)
#'  Then, for every cell-type, genes considered over-represented are saved in a list
#'
#'
#' @rdname heatmap_generation
#' @name heatmap_generation
#'
#' @param genesIn A list of gene symbols (all caps) to have their cell type enrichment
#' @param comp The name of the comparison
#' @param Cex  the size of the genes in the column label for the heatmap
#' @param rd_path The directory to R data files -- if they are not in this directory, then the files will be downloaded.
#' @param isMax If you are taking the single best CT marker (T/F) -- TRUE not reccomended
#' @param isPvalue If the signature matrix is raw p-value (T/F) -- TRUE not reccomended 
#' @param cellTypes Colnames of the cell-types you will extract (passed to extract_genes_cell)
#' @param pVal The level of association a gene has within a cell type (passed to extract_genes_cell)
#' @param isMax only take the top most DEG for each cell-type (T/F) -- TRUE not reccomended
#' @param isBackground If the heatmap is from the entire signature matrix or just the inputted gene list (T/F). isBackground == TRUE is used for internal.
#' @param refence Path to signature matrix or the signature matrix itself.
#' @param which_species Species of gene symbols -- "human" or "mouse" 
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#'
#' @return \code{heatmap_generation} A heatmap/barplot of p-value or odds-ratio of cell-type specific genes intersecting with the gene list. A list of genes that do/don't intersect with the signature matrix as well as a list of which cell-type these over-represented genes live in. \cr
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
#' # load in signature matrices
#' load("~/scMappR/data/Preoptic_region_example.rda")
#' POA_generes <- POA_example$POA_generes
#' POA_OR_signature <- POA_example$POA_OR_signature
#' POA_Rank_signature <- POA_example$POA_Rank_signature
#' # data(Preoptic_region_example)
#' Signature <- POA_Rank_signature
#'  rowname <- get_gene_symbol(Signature)
#'  rownames(Signature) <- rowname$rowname
#'  genes <- rownames(Signature)[1:100]
#'  heatmap_test <- heatmap_generation(genesIn = genes, "scMappR_test", reference = Signature, which_species = "mouse")
#'  
NULL
#' @rdname heatmap_generation
#' @export
#' 
heatmap_generation <- function(genesIn, comp, cex = 0.8, rd_path = "~/scMappR/data", cellTypes = "ALL", pVal = 0.01, isPval=TRUE, isMax =F,  isBackground = F,reference = "C:/Users/Dustin Sokolowski/Desktop/romanov_wilcoxon_test_2.RData",  which_species = "human", toSave = FALSE) {
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
  
  print(paste0("assumes inputted gene symbols are of ", which_species, " origin."))
  if(class(genesIn) != "character") stop("Error: inputted gene list is not character vector")
 
  if(class(reference) == "data.frame") reference <- as.matrix(reference)
  if(class(reference) == "matrix") { # then you have a custom signature matrix or it's pre-loaded into R
    wilcoxon_rank_mat_t <- reference  
  } else {
    if(class(reference) != "character") stop("Error: reference was not in matrix format or name of signature matrix RData file")
    load(reference) 
    
    wilcoxon_rank_mat_t <- wilcoxon_rank_mat_t[!duplicated(rownames(wilcoxon_rank_mat_t)),]
    if(length(grep("-", rownames(wilcoxon_rank_mat_t))) / length(rownames(wilcoxon_rank_mat_t)) > 0.75) {  
      print("Detected signature matrix from scMappR catelogue", quote = F)
      RN_2 <- get_gene_symbol(wilcoxon_rank_mat_t)
      
      rownames(wilcoxon_rank_mat_t) <- RN_2$rowname
      
    }
    if(which_species != RN_2$species ) { # convert background species to your species (internal only)
      
      warning(paste0("Species in signature matrix, ",RN_2$species, " is not the same as the inputted gene species ", which_species,"."))
      warning(paste0( "converting gene symbols in background from ",RN_2$species, " to ", which_species))
      
      
      
      thefiles <- list.files(path = rda_path, "bioMart_ortholog_human_mouse.rda")
      
      
      if(length(thefiles) == 0) {
        warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
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
        load(paste0(rda_path,"/bioMart_ortholog_human_mouse.rda"))
      }
      
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
    print("cell-types not uniquely named, appending tag to make all cell-types unique")
    colnames(wilcoxon_rank_mat_t) <- paste0(colnames(wilcoxon_rank_mat_t),"_tag",1:ncol(wilcoxon_rank_mat_t))
    
  }
  wilcoxon_rank_mat_t <- wilcoxon_rank_mat_t[!duplicated(rownames(wilcoxon_rank_mat_t)),] # just incase remove genes with the same name
  
  genesInter <- intersect(genesIn, rownames(wilcoxon_rank_mat_t)) # get genes with some DE
  whichGenesInter <- which(rownames(wilcoxon_rank_mat_t) %in% genesInter)
  genes_noInter <- genesIn[!(genesIn %in% genesInter)]
  geneHeat <- c()
  if(isBackground == TRUE) {
    print("preferential expression in background set: ")
  }
  if(isBackground == FALSE) {
    print("preferential expression in inputted gene list: ")
  }
  if(length(whichGenesInter) == 0) { # if no genes inputted are preferentially expressed
    print("No input genes are preferentially expressed")
    return("No_downstram_analysis")
  }
  if(length(whichGenesInter) == 1) { # if only 1 gene is
    print("One input gene is preferentially expressed")
    if(toSave == TRUE) {
    pdf(paste0(comp,"_barplot.pdf"))
    graphics::barplot(wilcoxon_rank_mat_t[whichGenesInter,], las = 2, main = rownames(wilcoxon_rank_mat_t)[whichGenesInter])
    dev.off()
    } else {
      warning("toSave = F and threfore plots are not allowed to be saved. I would reccomend allowing it to be true.")
    }
    return("No_downstream_analysis")
  }
  if(length(whichGenesInter) > 1) { # if > 1 genes are
    print("at least one input gene is preferentially expressed")
    # make the heatmap
    if(toSave == TRUE) {
    myheatcol <- grDevices::colorRampPalette(c("lightblue", "white", "orange"))(256)
    pdf(paste0(comp,"_heatmap.pdf"))
    gplots::heatmap.2(wilcoxon_rank_mat_t[whichGenesInter,], Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
    dev.off()
    geneHeat <- wilcoxon_rank_mat_t[whichGenesInter,]
    preferences <- extract_genes_cell(geneHeat, cellTypes = cellTypes, val = pVal, isMax = isMax, isPvalue = isPval)

    save(preferences, file = paste0(comp,"_preferences.RData"))
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
