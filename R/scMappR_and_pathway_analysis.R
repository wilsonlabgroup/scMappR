#' Generate STV, visualize, and enrich.
#' 
#' This function generates STVs, visualies them in a heatmap, and completes pathway enrichment of STVs and bulk gene list.
#'
#' This function generates STVs for every cell-type (see deconvolute_and_contextualize) as well as the relative cell-type proportions (which will be reutrned and pushed through).
#' Then, it generates heatmaps of all STVs, STVs overlapping with the signature matrix, the signature matrix, the signature matrix overlapping with STVs.
#' Then, if you have WIFI, it will complete g:ProfilR of the reordered STVs as well as a the ordered list of genes.
#' This function is a wrapper for deconvolute_and_contextualize, as well as pathway_enrich_internal.
#' 
#' @rdname scMappR_and_pathway_analysis
#' @name scMappR_and_pathway_analysis
#' 
#' @param count_file Normalized RNA-seq count matrix where rows are gene symbols and columns are individuals. Either the object tself of the path of a TSV file
#' @param signature_matrix Signature matrix (reccommended odds ratios) of cell-type specificity of genes. Either the object itself or a pathway to an RData file containing an object named "wilcoxon_rank_mat_or" -- generally internal.
#' @param DEG_list An object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted P-value, and the third the log2FC path to a tsv file containing this info is also acceptable.
#' @param case_grep Tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases.
#' @param control_grep tag in the column name for controls (i.e. samples representing downregulated OR an index of controls.
#' @param max_proportion_change Maximum cell-type proportion change -- may be useful if there are many rare cell-types.
#' @param print_plots Whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed. The same name of the plots will be completed for top pathways
#' @param plot_names The prefix of plot pdf files.
#' @param output_directory The name of the directory that will contain off of the analysis.
#' @param theSpecies -9 if using a precomputed count matrix from scMappR, human, mouse, or a specied directly compatible with g:Profiler. Removes ensembl symbols if appended
#' @param sig_matrix_size Number of genes in signature matrix for cell-type deconvolution.
#' @param drop_unkown_celltype Whether or not to remove "unknown" cell-types from the signature matrix.
#' @param WIFI Whether you have stable WIFI (T/F).
#' @param up_and_downregulated Whether you are additionally splitting up/downregulated genes (T/F).
#' @param gene_label_size The size of the gene label on the plot.
#' @param number_genes The number of genes to cut-off for pathway analysis (good with many DEGs)
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#' 
#' @return \code{scMappR_and_pathway_analysis} A directory with: STVs in RData file, Cell Type proportions (RData file), cell-type proportions leave one out (RData file), heatmap of STVs (all), heatmap of STVs (within signature), heatmap of signature (all), heatmap of signature (overlapping with DEG_list), Pathway enrichment for DEG list(all), RData file and Biological Processes, Pathway enrichment of STVs for each cell-type, RData file and biological processes. \cr
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
#' \notrun {
#' 
#' load("~/scMappR/data/PBMC_scMappR_and_pathway_analysis_example.rda")
#' # data(PBMC_scMappR_and_pathway_analysis)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in

#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, bulk_DE_cors, case_grep = case_grep, control_grep = control_grep, rda_path = "", max_proportion_change = 10, print_plots = T, plot_names = "tst1", theSpecies = "human", output_directory = "tester", sig_matrix_size = 3000, up_and_downregulated = F, internet = F)
#' 
#' }
#' 
NULL
#' @rdname scMappR_and_pathway_analysis
#' @export
#' 
scMappR_and_pathway_analysis <- function(  count_file,signature_matrix, DEG_list, case_grep, control_grep, rda_path = "", max_proportion_change = -9, print_plots=T, plot_names="scMappR",theSpecies = "human", output_directory = "scMappR_analysis",sig_matrix_size = 3000, drop_unkown_celltype = TRUE, internet = TRUE, up_and_downregulated = FALSE, gene_label_size = 0.4, number_genes = -9, toSave=FALSE) {
  
  
  
  # This function generates STVs for every cell-type (see deconvolute_and_contextualize) as well as the relative cell-type proportions (which will be reutrned and pushed through)
  # Then, it generates heatmaps of all STVs, STVs overlapping with the signature matrix, the signature matrix, the signature matrix overlapping with STVs
  # Then, if you have WIFI, it will complete g:ProfilR of the reordered STVs as well as a the ordered list of genes.
  # This function is a wrapper for deconvolute_and_contextualize, as well as pathway_enrich_internal  
  
  # Args:
  # count_file: Normalized RNA-seq count matrix where rows are gene symbols and columns are individuals
  # either the object tself of the path of a TSV file
  # signature_matrix: signature matrix (reccommended odds ratios) of cell-type specificity of genes
  # either the object itself or a pathway to an RData file containing an object named "wilcoxon_rank_mat_or" -- generally internal
  # DEG_list
  # an object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted P-value, and the third the log2FC
  # path to a tsv file containing this info is also acceptable
  #case_grep
  # tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases
  #control_grep
  # tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases
  #max_proportion_change
  # maximum cell-type proportion change -- may be useful if there are many rare cell-types
  # print_plots 
  # whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed
  # The same name of the plots will be completed for top pathways
  #plot_names
  # if plots are being printed, the prefix of their pdf files
  # output_directory
  # the name of the directory that will contain off of the analysis
  #theSpecies
  # -9 if using a precomputed count matrix from scMappR, human otherwise.
  # removes ensembl symbols
  #sig_matrix_size 
  # number of genes in signature matrix for cell-type deconvolution
  # drop_unkown_celltype
  # whether or not to remove "unknown" cell-types from the signature matrix
  # WIFI: Whether you have stable WIFI -- T/F
  # up_and_downregulated: whether you are splitting up/downregulated genes -- T/F
  # gene_label_size = the size of the gene label on the plot
  
  # Returns: 
  # a directory with:
  # STVs in RData file
  # Cell Type proportions (RData file)
  # cell-type proportions leave one out (RData file)
  # heatmap of STVs (all)
  # heatmap of STVs (within signature)
  # heatmap of signature (all)
  # heatmap of signature (overlapping with DEG_list)
  # Pathway enrichment for DEG list(all)
  # RData file and Biological Processes
  # Pathway enrichment of STVs for each cell-type
  # RData file and biological processes
  
  # load in the count matrix

  
  theSpecies <- tolower(theSpecies)
  if(class(count_file) == "character") {
    norm_counts_i <- read.table(count_file, header = T, as.is = T, sep = "\t")
  } else {
    norm_counts_i <- count_file
  }
  # get the background for later pathway analysis
  background_genes <-  rownames(norm_counts_i)  
  # list of differential expression
  if(class(DEG_list) == "character") {
    DEGs <- read.table(DEG_list, header = F, as.is = T, sep = "\t")
  } else {
    DEGs <- as.data.frame(DEG_list)
  }
  
  colnames(DEGs) <- c("gene_name", "padj", "log2fc")
  rownames(DEGs) <- DEGs$gene_name
  if(class(case_grep) != "character" | length(case_grep) > 1) {
    print("Assuming that case_grep and control_grep are indeces of 'case' and 'control'.", quote = F)
    print("Appending 'scMappR_case' to cases and 'scMappR_control to controls.", quote = F  )
    colnames(count_file)[case_grep] <- paste0("scMappR_case_", colnames(count_file)[case_grep])
    colnames(count_file)[control_grep] <- paste0("scMappR_control_", colnames(count_file)[control_grep])
    control_grep <- "scMappR_control"
    case_grep <- "scMappR_case"
  }
  if(class(signature_matrix) == "character") {
    # assuming that the signature matrix is an RData file containg an object named wilcoxon_rank_mat_or (internal, generally)
    load(signature_matrix)
    signature_matrix <- wilcoxon_rank_mat_or
  }
  if(length(grep("-", rownames(signature_matrix))) / length(rownames(signature_matrix)) > 0.75 ) {  
    warning("More than 50 genes contian '-', and the signature matrix is considered internal")
    load(paste0(rda_path,"/bioMart_ortholog_human_mouse.rda"))
    # data(bioMart_ortholog_human_mouse)    
    RN_2 <- get_gene_symbol(signature_matrix)
    
    rownames(signature_matrix) <- RN_2$rowname
    print(head(signature_matrix))
    #stop("testing")
    internal_species <- RN_2$species
    if(internal_species != theSpecies) {
      warning(paste0("the species from the signature matrix, ", internal_species, ", does not equal the initially inputted species, ", theSpecies, ". Converting gene symbols of 1-1 orthologs"))
      RN <- rownames(signature_matrix) # gene symbols

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
      
      # data(bioMart_ortholog_human_mouse) 
      
      rownames(bioMart_orthologs) <- bioMart_orthologs[,internal_species]
      intersected_genes <- intersect(RN, rownames(bioMart_orthologs)) # gene sybmols in signature
      signature_matrix1 <- signature_matrix[intersected_genes,]
      bioMart_orthologs1 <- bioMart_orthologs[intersected_genes,]
      rownames(signature_matrix1) <- bioMart_orthologs1[,theSpecies] # replacing rownames with the species you want
      print(dim(signature_matrix1))
      print(head(rownames(signature_matrix1)))
      signature_matrix <- signature_matrix1
    }
  }
  
  # Identify scMappR Transformed Values 
  STVs <- deconvolute_and_contextualize(count_file, signature_matrix, DEG_list, case_grep , control_grep, max_proportion_change = max_proportion_change, print_plots = print_plots, plot_names = plot_names, theSpecies = theSpecies, sig_matrix_size = sig_matrix_size, drop_unkown_celltype = drop_unkown_celltype)
  
  print(paste0("Making scMappR output directory named", output_directory, "."))
  if(toSave == FALSE) {
    warning("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended. Returning STVs and no pathway analysis.")
    print("toSave == FALSE therefore files cannot be saved. Switching toSave = TRUE is strongly reccomended. Returning STVs and no pathway analysis.", quote = F)
    return(STVs)    
  }
  dir.create(output_directory)
  
  scMappR_vals <- STVs$scMappR_transformed_values # scMappR values
  print(scMappR_vals)
  save(scMappR_vals, file = paste0(output_directory, "/",plot_names, "_STVs.RData"))
  
  cell_proportions_all <- STVs$cellType_Proportions # all gene CT proportion
  save(cell_proportions_all, file = paste0(output_directory, "/",plot_names, "_celltype_proportions.RData"))
  
  leave_one_out_proportions <- STVs$leave_one_out_proportions # leave one out avg CT proportions
  save(cell_proportions_all, file = paste0(output_directory, "/",plot_names, "_leaveOneOut_gene_proportions.RData"))
  
  signature_mat <- STVs$processed_signature_matrix # processed_signaure_matrix
  save(signature_mat, file = paste0(output_directory, "/",plot_names, "_leaveOneOut_gene_proportions.RData"))
  if(nrow(DEG_list) == 1) {
    warning("You only have 1 DEG, no heatmaps can be made. Returning STV")
    print("You only have 1 DEG, no heatmaps can be made. Returning STV",quote = F)
    return(scMappR_vals)
    }
  myheatcol <- grDevices::colorRampPalette(c("lightblue", "white", "orange"))(256)
  
  # generate heatmaps for DEGs
  cex = gene_label_size
  pdf(paste0(output_directory, "/", plot_names,"_STVs_all_genes_heatmap.pdf"))
  gplots::heatmap.2(as.matrix(scMappR_vals), Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
  dev.off()
  
  pdf(paste0(output_directory, "/",plot_names,"_signature_matrix_all_genes_heatmap.pdf"))
  gplots::heatmap.2(as.matrix(signature_mat), Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
  dev.off()
  
  celltype_preferred_degs <- intersect(rownames(scMappR_vals), rownames(signature_matrix)) # intersect DEGs and Genes in signature matrix
  
  if(length(celltype_preferred_degs) < 3) {
    warning("Fewer than 3 genes are both De and in the signature matrix. Therefore, these heatmaps will not be generated. Furthermore, there is insufficient re-ranking of genes for different pathway analyses to be neccesary. Therefore, here, STVs are more representative of a scaling factor for each cell-type.")
    print("Fewer than 3 genes are both De and in the signature matrix. Therefore, these heatmaps will not be generated. Furthermore, there is insufficient re-ranking of genes for different pathway analyses to be neccesary. Therefore, here, STVs are more representative of a scaling factor for each cell-type.", quote = F)
    
    return(list(STVs = STVs))

  } else {
    
    # generating the heatmaps for STVs and signature matrix odds ratios that overlap with one another
    
    pdf(paste0(output_directory, "/", plot_names,"_STVs_signature_genes_heatmap.pdf"))
    gplots::heatmap.2(as.matrix(scMappR_vals[celltype_preferred_degs,]), Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
    dev.off()
    
    pdf(paste0(output_directory, "/",plot_names,"_signature_matrix_DEGs_heatmap.pdf"))
    gplots::heatmap.2(as.matrix(signature_mat[celltype_preferred_degs,]), Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
    dev.off()
    
  }
  if(internet == FALSE) {
    warning("There is not a reported stable internet (WIFI = FALSE) and therefore pathway analysis with g:Prof")
    return("Done!")
  }
  up_and_down_together <- pathway_enrich_internal(  DEGs, theSpecies, scMappR_vals, background_genes, output_directory, plot_names, number_genes, toSave=TRUE)
  if(up_and_downregulated == TRUE)  {
    print("Splitting genes by up- and down-regulated and then repeating analysis", quote = F)
    rownames(DEGs) <- DEGs$gene_name
    upGenes <- DEGs$gene_name[DEGs$log2fc > 0]
    downGenes <- DEGs$gene_name[DEGs$log2fc < 0]
    upDir <- paste0(output_directory,"/upregulated")
    downDir <- paste0(output_directory, "/downregulated")
    dir.create(upDir)
    dir.create(downDir)
    
    upDEGs <- DEGs[upGenes,]
    downDEGs <- DEGs[downGenes,]
    upSTVs <- scMappR_vals[upGenes,]
    DownSTVs <- scMappR_vals[downGenes,]
    print("Pathway analysis of upregulated genes")
    up_only <- pathway_enrich_internal(  upDEGs, theSpecies, upSTVs, background_genes, upDir, plot_names, number_genes, toSave=TRUE)
    print("Pathway analysis of downregulated genes")
    down_only <- pathway_enrich_internal(  downDEGs, theSpecies, DownSTVs, background_genes, downDir, plot_names, number_genes, toSave=TRUE)    
  }
  
  return(list(STVs = STVs, paths = up_and_down_together$biological_pathways, TFs = up_and_down_together$transcription_factors))
}
