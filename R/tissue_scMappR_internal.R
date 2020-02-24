#' Gene List Visualization and Enrichment (Internal)
#'
#' This function loops through every signature matrix in a particular tissue and generates heatmaps, cell-type preferences, and co-enrichment.
#' 
#' This function takes a list of genes and a tissue that is contained in current signature matrices before and generating heatmaps of cell-type preferences.
#' It then completes cell-type enrichment of each individual cell-type, then, if more than two cell-types are signficiantly enriched, co-enrichemnt 
#' of those enriched cell-types is then computed.
#'
#'
#' @rdname tissue_scMappR_internal
#' @name tissue_scMappR_internal
#'
#' @param gene_list A list of gene symbols, mouse or human.
#' @param species "mouse", "human" or "-9" if using a precomputed signature matrix.
#' @param tissue  Name of the tissue in "get_tissues".
#' @param cluster 'Pval' or 'OR' depending on if you wantto cluster odds ratios or pvalues of CT preferences.
#' @param genecex The size of the gene names of the rows in the heatmap.
#' @param raw_pval  If the inputted signature matrix are raw (untransformed) p-values -- recommended to generate rank first (T/F).
#' @param rda_path Path to the .rda file containing all of the signature matrices.
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' @param output_directory If toSave = TRUE, the name of the output directory that would be built.
#' @param drop_unkown_celltype Whether or not to remove "unknown" cell-types from the signature matrix (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#' 
#' @return List with the following elements:
#' \item{background_heatmap}{Data frame of the entire gene by cell-type signature matrix inputted.}
#' \item{gene_list_heatmap}{Data frame of inputted signature matrix subsetted by input genes.}
#' \item{single_celltype_preferences}{Data frame of enriched cell-types.}
#' \item{group_celtype_preference}{Data frame of groups of cell-types enriched by the same genes.}
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
#' 
#' data(POA_example) # region to preoptic area
#' Signature <- POA_example$POA_Rank_signature # signature matrix
#' rowname <- get_gene_symbol(Signature) # get signature
#' rownames(Signature) <- rowname$rowname
#' genes <- rownames(Signature)[1:60]
#' rda_path1 = "" # data directory (if it exists)
#' 
#' # set toSave = TRUE and path = output directory of your choice
#' internal <- tissue_scMappR_internal(genes, "mouse", output_directory = "scMappR_TesInternal",
#'                                    tissue = "hypothalamus", toSave = FALSE) 
#' 
#'  }
#' @export
#' 
tissue_scMappR_internal <- function(gene_list,species, output_directory, tissue,rda_path = "", cluster = "Pval", genecex = 0.01, raw_pval = FALSE, path = NULL, toSave = FALSE, drop_unkown_celltype = FALSE) {
  
  # This function takes a list of genes and a tissue available in 'get_tissues' function and generates heatmaps of cell-type preferences
  # it then completes cell-type enrichment of each individual cell-type, then, if more than two cell-types are signficiantly enriched, co-enrichemnt 
  # of those cell-types
  # Args:
  # gene_list: a list of gene symbols, mouseor human
  # species: "mouse", "human" or "-9" if using a precomputed signature matrix
  # tissue: name of the tissue in "get_tissues"
  # cluster: 'Pval' or 'OR' depending on if you wantto cluster odds ratios or pvalues of CT preferences.
  # genecex: the size of the genes in the rows of the heatmap
  # raw_pval: if the inputed signature matrix are raw (untransformed) Pvalues -- reccomended to generate rank first
  # Returns: generates a directory containing
  # heatmap for signature matrix and matrix intersecting with heatmap
  # RData file containing cell-type preferences
  # tsv files with gene set enrichment for single cell-types and co-enrichment
  
  if(!is.character(gene_list)) {
    stop("The 'gene_list' argumet must be a character vector of gene symbols.")
  }
  if(!(species %in% c("mouse", "human"))) {
    stop("the 'species' argument must be 'mouse' or 'human' (case sensitive). ")
  }
  

  if(!is.character(output_directory)) {
    stop("output_directory must be of class character.")
  }
  
  if(!is.numeric(genecex)) {
    stop("genecex must be of class numeric.")
  }
  if(genecex < 0) {
    warning("making genecex low positive value as it is currently < 0.")
    genecex <- 0.001
  }

  
  if(!is.character(cluster)) {
    stop("cluster must be of class character.")
  }
  if(!is.character(rda_path)) {
    stop("rda_path must be of class character.")
  }

  if(all(is.logical(toSave), is.logical(raw_pval), is.logical(drop_unkown_celltype)) == FALSE) {
    stop("toSave and raw_pval and drop_unknown_celltype must be of class logical.")
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  RankValueSignature <- "" # empty for cran
  OddsRatioSignature <- "" # 
  outDir <- output_directory
  if(cluster == "Pval") {
    clust <- "p_val"
    
  } else {
    clust <- "odd_ratio"
  }
  if(raw_pval == TRUE) {
    warning("Generating heatmap of Raw P-values, the lower the heat the more significant. Would reccomend switching to ranks insteads (-log10(Pval) * sign(Fold-Change))")
  }
  if(clust == "p_val") {
    
    thefiles <- list.files(path = rda_path, "Signature_matrices_Pval.rda")
    
    
    if(length(thefiles) == 0) {
      warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
      #
      datafile <- "Signature_matrices_Pval.rda"
      metafile <- paste0(datafile)
      url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      #
    } else {
      load(paste0(rda_path,"/Signature_matrices_Pval.rda"))
    }
    
    scMappR_list <- RankValueSignature
  } else {
    
    thefiles <- list.files(path = rda_path, "Signature_matrices_OR.rda")
    
    
    if(length(thefiles) == 0) {
      
      warning(paste0("Cell-marker databases are not present in ", rda_path, " downloading and loading data."))
      #
      datafile <- "Signature_matrices_OR.rda"
      metafile <- paste0(datafile)
      url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      #
    } else {
      load(paste0(rda_path,"/Signature_matrices_OR.rda"))
    }
    
    
    
    scMappR_list <- OddsRatioSignature
  }
  hm <- names(scMappR_list) # get all the available pvalue signautres
  hm_tissue <- grep(toupper(tissue), toupper(hm))
  if(length(hm_tissue) == 0) {
    stop(paste0(tissue, " is not currently an available tissue for scMappR. Please check spelling or try a different tissue."))
  }
  input_studies <- hm[hm_tissue]
  study_names <- input_studies
  if(toSave == TRUE) {
    dir.create(paste0(path,"/",outDir))
  } else {
    warning("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.")
    message("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.", quote = FALSE)
    
  }
  single_cell_studies <- list()
  for(i in 1:length(study_names)) {
    
    wilcoxon_rank_mat_t <- scMappR_list[[study_names[i]]]
    
    if(length(grep("-", rownames(wilcoxon_rank_mat_t))) / length(rownames(wilcoxon_rank_mat_t)) > 0.75) {  
      message("Detected signature matrix from scMappR catelogue", quote = FALSE)
      RN_2 <- get_gene_symbol(wilcoxon_rank_mat_t)
      
      rownames(wilcoxon_rank_mat_t) <- RN_2$rowname
      
    }
    if(species != RN_2$species ) { # convert background species to your species (internal only)
      
      warning(paste0("Species in signature matrix, ",RN_2$species, " is not the same as the inputted gene species ", species,"."))
      warning(paste0( "converting gene symbols in background from ",RN_2$species, " to ", species))
      
      thefiles <- list.files(path = rda_path, "ortholog_human_mouse.rda")
      
      if(length(thefiles) == 0) {
        warning(paste0("Cell-marker preferences are not present in ", rda_path, " downloading and loading data."))
        #
        datafile <- "bioMart_ortholog_human_mouse.rda"
        metafile <- paste0(datafile)
        url <- paste0("https://github.com/DustinSokolowski/scMappR_Data/blob/master/", 
                      metafile, "?raw=true")
        destfile <- file.path(tempdir(), metafile)
        downloader::download(url, destfile = destfile, mode = "wb")
        load(destfile)
      } else {
        
        load(paste0(rda_path,"/bioMart_ortholog_human_mouse.rda"))
      }
      
      
      
      rownames(bioMart_orthologs) <- bioMart_orthologs[,RN_2$species]
      RN <- rownames(wilcoxon_rank_mat_t)
      
      intersected_genes <- intersect(RN, rownames(bioMart_orthologs)) # gene sybmols in signature
      wilcoxon_gene1 <- wilcoxon_rank_mat_t[intersected_genes,]
      
      bioMart_orthologs1 <- bioMart_orthologs[intersected_genes,]
      rownames(wilcoxon_gene1) <- bioMart_orthologs1[,species] # replacing rownames with the species you want
      wilcoxon_rank_mat_t <- wilcoxon_gene1
      
      background_genes <- rownames(wilcoxon_rank_mat_t)
      theSpecies <- species
    } else {
      sym <- RN_2 # get the species and shave on ensembl symbol
      
      background_genes <- sym$rowname # genes of background
      theSpecies <- sym$species # name of species
    }
    study_ref <- wilcoxon_rank_mat_t
    
    unknown <- grep(toupper("unknown"),toupper(colnames(study_ref)))
    if(length(unknown) > 0 & drop_unkown_celltype == TRUE) {
      message("Removing unknown cell-types")
      study_ref <- study_ref[,-unknown]
    } else {
      study_ref <- study_ref
    }
    
    background_genes <- rownames(study_ref)

    background_heatmap <- heatmap_generation(background_genes, comp = paste0(outDir, "/", study_names[i],"_background"), reference = study_ref, isBackground = TRUE, cex = genecex, which_species = species, isPval = raw_pval, toSave=toSave, path = path)  
    # get the heatmap of all of the genes in the signature matrix
    message("Number of DEGs that are cell-type markers in current signature matrix: ", quote = FALSE)
    message(length(intersect(gene_list, rownames(study_ref))))
    theL <- length(intersect(gene_list, rownames(study_ref)))
    if(theL < 3) {
      message(paste0("Your gene list contains fewer than 3 overlapping genes with ",study_names[i],". Therefore no heatmap was saved and enrichment cannot be done."), quote = FALSE)
      message(paste0("Subsetted CT marker preferences of these genes are saved in ",paste0(outDir, "/", study_names[i],"_genelist")), quote = FALSE)
      message(intersect(gene_list, rownames(study_ref)))
      subsetted_genes <- study_ref[intersect(gene_list, rownames(study_ref)),]
      if(toSave==TRUE) {
      save(subsetted_genes, file = paste0(path, "/",outDir, "/", study_names[i],"_subsetted.RData"))
      } else {
        warning("Cannot save preferences of subsetted genes as toSave = FALSE")
      }
      next
     }

    gene_list_heatmap <- heatmap_generation(gene_list, comp = paste0(outDir, "/", study_names[i],"_genelist"), reference = study_ref, cex = genecex, which_species = species, isPval = raw_pval, toSave = toSave, path = path)
    
    # get the heatmap of genes overlapping with the signature matrix and the inputted gene list
    if(is.character(gene_list_heatmap)) {
      message("not enough genes were present to do downsteam analysis in: ")
      message(study_names[i])
      single_cell_studies[[i]] <- gene_list_heatmap
      next
    }
    singleCTpreferences <- single_gene_preferences(gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory, toSave = toSave, path = path)
    # interrogate the enrichment of every cell-type within the signature matrix
    sig <- singleCTpreferences[singleCTpreferences$pFDR < 0.05,]
    sig$Odds_Ratio <- toNum(sig$Odds_Ratio)
    sig <- sig[sig$Odds_Ratio > 1,]
    
    if(nrow(sig) <2 ) { 
      #if there are fewer than 2 enriched cell-types
      message("co-enrichment cannot be measured as one or fewer CTs are enriched")
      coCTpreferences <- "co-enrichment cannot be measured as one or fewer CTs are enriched"
      output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
      single_cell_studies[[i]] <- output 
      next
    }
    sig <- sig[order(sig$pFDR),]
    coCTpreferences <- coEnrich(sig, gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory, toSave = toSave, path = path)
    # complete co-enrichment of up to 5 cell-types
    output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
    
    single_cell_studies[[i]] <- output 
  }
  names(single_cell_studies) <- study_names
  return(single_cell_studies)
}

