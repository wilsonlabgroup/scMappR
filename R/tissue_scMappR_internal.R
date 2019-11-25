#' Gene List Visualization and Enrichment (Internal)
#'
#' This function loops through every signature matrix in a particular tissue and generates heatmaps, cell-type preferences, and co-enrichment.
#' 
#' This function takes a list of genes and a tissue available in 'get_tissues' function and generates heatmaps of cell-type preferences
#' it then completes cell-type enrichment of each individual cell-type, then, if more than two cell-types are signficiantly enriched, co-enrichemnt 
#' of those cell-types.
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
#' @param raw_pval  If the inputed signature matrix are raw (untransformed) Pvalues -- reccomended to generate rank first.
#' @param rda_path Path to the rda file containing all of the signature matrices.
#' @param toSave Allow scMappR to write files in the current directory (T/F)
#'
#' @return \code{tissue_scMappR_internal} A list containing the entire signature matrix, the matrix subsetted for your genes, enrichment of each cell-type, and co-enrichment. \cr
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
#' # load in signature matrices
#' load("~/scMappR/data/Preoptic_region_example.rda")
#'  POA_generes <- POA_example$POA_generes
#'  POA_OR_signature <- POA_example$POA_OR_signature
#'  POA_Rank_signature <- POA_example$POA_Rank_signature
#' 
#' # data(Preoptic_region_example)
#' Signature <- POA_Rank_signature
#'  rowname <- get_gene_symbol(Signature)
#'  rownames(Signature) <- rowname$rowname
#'  genes <- rownames(Signature)[1:200]
#'  # Assuming Signature_matrices_pVal.rda is in you "~/Documents/scMappR/data" directory
#'  # rda_path1 = "~/Documents/scMappR/data"
#'  internal <- tissue_scMappR_internal(genes,"mouse",output_directory = "scMappR_Test", tissue = "hypothalamus",rda_path = rda_path1)
#'  }
NULL
#' @rdname tissue_scMappR_internal
#' @export
#' 
tissue_scMappR_internal <- function(gene_list,species, output_directory, tissue,rda_path, cluster = "Pval", genecex = 0.01, raw_pval = FALSE, toSave = FALSE) {
  
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
  input_studies <- hm[hm_tissue]
  study_names <- input_studies
  if(toSave == TRUE) {
    dir.create(outDir)
  } else {
    warning("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.")
    print("toSave == FALSE and therefore a directory cannot be made. Switching toSave = TRUE is reccomended.", quote = F)
    
  }
  single_cell_studies <- list()
  for(i in 1:length(study_names)) {
    
    wilcoxon_rank_mat_t <- scMappR_list[[study_names[i]]]
    
    if(length(grep("-", rownames(wilcoxon_rank_mat_t))) / length(rownames(wilcoxon_rank_mat_t)) > 0.75) {  
      print("Detected signature matrix from scMappR catelogue", quote = F)
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
      sym <- get_gene_symbol(wilcoxon_rank_mat_t) # get the species and shave on ensembl symbol
      
      background_genes <- sym$rowname # genes of background
      theSpecies <- sym$species # name of species
    }
    study_ref <- wilcoxon_rank_mat_t
    background_genes <- rownames(study_ref)
    background_heatmap <- heatmap_generation(background_genes, comp = paste0(outDir, "/", study_names[i],"_background"), reference = study_ref, isBackground = TRUE, cex = genecex, which_species = species, isPval = raw_pval, toSave=toSave)  
    # get the heatmap of all of the genes in the signature matrix
    print("Number of DEGs that are cell-type markers in current signature matrix: ", quote =F)
    print(length(intersect(gene_list, rownames(study_ref))))
    theL <- length(intersect(gene_list, rownames(study_ref)))
    if(theL < 3) {
      print(paste0("Your gene list contains fewer than 3 overlapping genes with ",study_names[i],". Therefore no heatmap was saved and enrichment cannot be done."), quote = F)
      print(paste0("Subsetted CT marker preferences of these genes are saved in ",paste0(outDir, "/", study_names[i],"_genelist")), quote = F)
      print(intersect(gene_list, rownames(study_ref)))
      subsetted_genes <- study_ref[intersect(gene_list, rownames(study_ref)),]
      if(toSave==TRUE) {
      save(subsetted_genes, file = paste0(outDir, "/", study_names[i],"_subsetted.RData"))
      } else {
        warning("Cannot save preferences of subsetted genes as toSave = FALSE")
      }
      next
     }
      
    gene_list_heatmap <- heatmap_generation(gene_list, comp = paste0(outDir, "/", study_names[i],"_genelist"), reference = study_ref, cex = genecex, which_species = species, isPval = raw_pval, toSave = toSave)
    
    # get the heatmap of genes overlapping with the signature matrix and the inputted gene list
    if(class(gene_list_heatmap) == "character") {
      print("not enough genes were present to do downsteam analysis in: ")
      print(study_names[i])
      single_cell_studies[[i]] <- gene_list_heatmap
      next
    }
    singleCTpreferences <- single_gene_preferences(gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory, toSave = toSave)
    # interrogate the enrichment of every cell-type within the signature matrix
    sig <- singleCTpreferences[singleCTpreferences$pFDR < 0.05,]
    sig$Odds_Ratio <- toNum(sig$Odds_Ratio)
    sig <- sig[sig$Odds_Ratio > 1,]
    
    if(nrow(sig) <2 ) { 
      #if there are fewer than 2 enriched cell-types
      print("co-enrichment cannot be measured as one or fewer CTs are enriched")
      coCTpreferences <- "co-enrichment cannot be measured as one or fewer CTs are enriched"
      output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
      single_cell_studies[[i]] <- output 
      next
    }
    sig <- sig[order(sig$pFDR),]
    coCTpreferences <- coEnrich(sig, gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory, toSave = toSave)
    # complete co-enrichment of up to 5 cell-types
    output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
    
    single_cell_studies[[i]] <- output 
  }
  names(single_cell_studies) <- study_names
  return(single_cell_studies)
}

