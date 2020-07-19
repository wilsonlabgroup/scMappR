#' two_method_pathway_enrichment
#'
#' Pathway analysis of each cell-type based on cell-type specificity and rank improvement by scMappR.
#'
#' This function re-ranks cwFoldChanges based on their absolute ct specificity scores (per-celltype) as well as their rank increase in cell-type specificity before completing an ordered pathway analysis. In the second method, only genes with a rank increase in cell-type specificity were included
#'  
#'
#' @rdname two_method_pathway_enrichment
#' @name two_method_pathway_enrichment
#'
#' @param DEGs Differentially expressed genes (gene_name, padj, log2fc).
#' @param theSpecies Human, mouse, or a charcter that is compatible with gProfileR.
#' @param scMappR_vals cell weighted Fold-changes of differentially expressed genes.
#' @param background_genes A list of background genes to test against. NULL assumes all genes in gprofiler gene set databases.
#' @param output_directory Path to the directory where files will be saved.
#' @param plot_names Names of output.
#' @param number_genes Number of genes to if there are many, many DEGs.
#' @param newGprofiler Whether to use gProfileR or gprofiler2 (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#' 
#' 
#' @return List with the following elements:
#' \item{rank_increase}{A list containing the degree of rank change between bulk DE genes and cwFold-changes. Pathway enrichment and tf enrichment of these reranked genes.}
#' \item{non_rank_increase}{list of DFs containing the pathway and TF enrichment of cwFold-changes.}
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
#' 
#' # load data for scMappR
#' data(PBMC_example)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' 
#' # calculate cwFold-changes
#' toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
#'                                       bulk_DE_cors, case_grep = case_grep,
#'                                       control_grep = control_grep, rda_path = "", 
#'                                       max_proportion_change = 10, print_plots = TRUE, 
#'                                       plot_names = "tst1", theSpecies = "human", 
#'                                       output_directory = "tester",
#'                                       sig_matrix_size = 3000, up_and_downregulated = FALSE, 
#'                                       internet = FALSE)
#' 
#' # complete pathway enrichment using both methods
#' twoOutFiles <- two_method_pathway_enrichment(bulk_DE_cors, "human",
#' scMappR_vals = toOut$cellWeighted_Foldchange, background_genes = rownames(bulk_normalized), 
#' output_directory = "newfun_test",plot_names = "nonreranked_", toSave = FALSE)
#' 
#' }
#'  
#' @export
#' 

two_method_pathway_enrichment <- function(DEGs, theSpecies, scMappR_vals, background_genes = NULL, output_directory = "output", plot_names = "reweighted", number_genes = -9,  newGprofiler = FALSE, toSave = FALSE, path = NULL) {

  scMappR_vals_class <- class(scMappR_vals)[1] %in% c("data.frame", "matrix")
  if(scMappR_vals_class[1] == FALSE) {
    stop("scMappR_vals must be a data.frame or matrix.")
  }
  bg_test <- all(!is.character(background_genes), !is.null(background_genes))[1]
  if(bg_test) {
    stop("background_genes must be a character vector of gene names.")
  }
  
  if(!is.character(theSpecies)) {
    stop("the species must be a character, human, mouse, or a species compatible with gprofiler.")
  }
  
  if(!is.character(output_directory)) {
    stop("output_directory must be a character.")
  }
  
  if(!is.character(plot_names)) {
    stop("plot_names must be a character, human, mouse, or a species compatible with gprofiler.")
  }
  if(!is.numeric(number_genes)) {
    message(number_genes)
    message(class(number_genes))
    stop("number_genes must be of class numeric.")
  }
  if(!is.logical(newGprofiler)) {
    stop("newGprofiler must be logical TRUE/FALSE.")
  }
  
  
  if(!toSave) {
    warning("toSave = FALSE and therefore scMappR is not allowed to print files and pictures. For this function to work fully, please set toSave = TRUE")
  }
  
  if(toSave) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  if(all(toSave, !(is.null(path)))[1]) {
    dir.create(paste0(path,"/",output_directory))
    dir.create(paste0(path,"/",output_directory,"/bp_reranked"))
    dir.create(paste0(path,"/",output_directory,"/tf_reranked"))
    dir.create(paste0(path,"/",output_directory,"/scatter_reranked"))
    
    
  }
  
  ##

  
  
  ##
  ### Analysis of the rank-increase of DE genes that increased in rank due to scMappR
  ##
  message("Completing traditional pathway analysis")
  ### Traditional pathway analysis of each cell-type -- completing CT specific pathway enrichment and not looking at how much more significant a gene was
  non_rank_increase <- scMappR::pathway_enrich_internal(DEGs, theSpecies, scMappR_vals, background_genes, output_directory, plot_names, number_genes,  newGprofiler, toSave, path)
  ##
  ##
  
  
  if(theSpecies == "human") species_bulk <- "hsapiens"
  if(theSpecies == "mouse") species_bulk <- "mmusculus"
  # Interrogating the pathway enrichment of the rank-increase in cell-type specificity
  celltype_path <- list() 
  for(i in 1:ncol(scMappR_vals)) {
    
    message(paste0("generating rank improvement DF for ", names(scMappR_vals)[i]))
    
    # generating the DF containing the rank order difference between cwFold-changes and DE genes
    scMappR_vals_CT <- scMappR_vals[order(abs(scMappR_vals[,i]), decreasing = TRUE),]

    rank_CT <- 1:length(rownames(scMappR_vals_CT))
    rank_DEGs <- 1:length(DEGs$gene_name)
    names(rank_CT) <- rownames(scMappR_vals_CT)
    names(rank_DEGs) <- DEGs$gene_name
    
    rank_CT <- rank_CT[intersect(names(rank_CT),names(rank_DEGs) )]
    rank_DEGs <- rank_DEGs[intersect(names(rank_CT),names(rank_DEGs) )]
    rank_change <- rank_CT - rank_DEGs
    ranks <- as.data.frame(cbind(rank_CT, rank_DEGs, rank_change))
    ranks <- ranks[order(ranks$rank_change),]  
    
    
    ranks$col <- "grey"
    ranks$col[ranks$rank_change < 0] <- "red"
    
    # Identifying gene list that improves in the context of ct specific pathway enrichment
    ranks_improve <- ranks[ranks$rank_change < 0,]

    # Scatterplot of DE genes
    saveWorks <- all(toSave, !is.null(path))[1]
    if(saveWorks) {    
    message(paste0("generating rank improvement scatterplot for ", names(scMappR_vals)[i]))
      
    grDevices::pdf( file = paste0(path,"/",output_directory,"/scatter_reranked/",colnames(scMappR_vals)[i],"_rank_scatter.pdf"))
    
    ranks$rank_DEGs <- as.numeric(ranks$rank_DEGs)
    ranks$rank_CT <- as.numeric(ranks$rank_CT)
    p <- ggplot2::ggplot(ranks, ggplot2::aes(x= rank_DEGs, y= rank_CT, colour= col)) + ggplot2::geom_point(ggplot2::aes(colour = rank_change), show.legend = FALSE) +   ggplot2::scale_colour_gradient2(mid = "grey", low = "darkgreen", high = "purple") + ggplot2::ggtitle(names(scMappR_vals)[i]) + ggplot2::theme_bw(base_size = 13) + ggplot2::xlab("Rank of Bulk DE genes") + ggplot2::ylab("Rank of cwFold-changes") + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), plot.title = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 12))
 
    print(p)
    
    grDevices::dev.off()
    
    # PNG
    grDevices::png( filename  = paste0(path,"/",output_directory,"/scatter_reranked/",colnames(scMappR_vals)[i],"_rank_scatter.png"))
    
    ranks$rank_DEGs <- as.numeric(ranks$rank_DEGs)
    ranks$rank_CT <- as.numeric(ranks$rank_CT)
    p <- ggplot2::ggplot(ranks, ggplot2::aes(x= rank_DEGs, y= rank_CT, colour= col)) + ggplot2::geom_point(ggplot2::aes(colour = rank_change), show.legend = FALSE) +   ggplot2::scale_colour_gradient2(mid = "grey", low = "darkgreen", high = "purple") + ggplot2::ggtitle(names(scMappR_vals)[i]) + ggplot2::theme_bw(base_size = 13) + ggplot2::xlab("Rank of Bulk DE genes") + ggplot2::ylab("Rank of cwFold-changes") + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), plot.title = ggplot2::element_text(size = 14), axis.title.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 12))
    
    print(p)
    
    grDevices::dev.off()
    
    
    }

    if(!newGprofiler) {
      
      if(is.null(background_genes)) background_genes <- ""
      
      #Old gprofiler, biological pathways
    scMappR_vals_specific1 <- gProfileR::gprofiler(rownames(ranks_improve), organism = species_bulk, ordered_query = TRUE, min_isect_size = 3, min_set_size = 15, max_set_size = 2000,src_filter = c("GO:BP", "REAC", "KEGG"), hier_filtering = "moderate", correction_method = "fdr")
    scMappR_vals_specific1$term_name <- scMappR_vals_specific1$term.name
    scMappR_vals_specific1$p_value <- scMappR_vals_specific1$p.value
    ordered_back_all <- scMappR_vals_specific1
    
      #Old gprofiler, transcription factors
    scMappR_vals_specific1 <- gProfileR::gprofiler(rownames(ranks_improve), organism = species_bulk, ordered_query = TRUE, min_isect_size = 3, min_set_size = 15, max_set_size = 5000,src_filter = c("TF"), hier_filtering = "moderate", correction_method = "fdr")
    #scMappR_vals_specific1 <- scMappR_vals_gprof[!(scMappR_vals_gprof$term.name %in% DEGs_gprof$term.name),]
    scMappR_vals_specific1$term_name <- scMappR_vals_specific1$term.name
    scMappR_vals_specific1$p_value <- scMappR_vals_specific1$p.value
    ordered_back_all_tf <- scMappR_vals_specific1
    
    
    } else {
      ordered_back_all <- gprofiler2::gost(query = rownames(ranks_improve), organism = species_bulk, ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "fdr",  custom_bg =background_genes, numeric_ns = "", sources = c("GO:BP", "KEGG", "REAC"))  
      
      if(is.null(ordered_back_all)) { # grpofiler2 returns null if nothing is significant, making compatile with plotBP and make_TF_barplot
        ordered_back_all <- as.data.frame(matrix(c(1," "),nrow = 1))
        colnames(ordered_back_all) <- c("p_value", "term_name")
        ordered_back_all <- ordered_back_all[-1,]
      } else {
        ordered_back_all <- ordered_back_all$result
        ordered_back_all <- ordered_back_all[ordered_back_all$term_size > 15 & ordered_back_all$term_size < 2000 & ordered_back_all$intersection_size > 2,]
      }  
      
      ordered_back_all_tf <- gprofiler2::gost(query = rownames(ranks_improve), organism = species_bulk, ordered_query = TRUE, significant = TRUE, exclude_iea = FALSE, multi_query = FALSE, measure_underrepresentation = FALSE, evcodes = FALSE, user_threshold = 0.05, correction_method = "fdr", custom_bg =background_genes, numeric_ns = "", sources = c("TF"))  
      
      if(is.null(ordered_back_all_tf)) { # grpofiler2 returns null if nothing is significant, making compatile with plotBP and make_TF_barplot
        ordered_back_all_tf <- as.data.frame(matrix(c(1," "),nrow = 1))
        colnames(ordered_back_all_tf) <- c("p_value", "term_name")
        ordered_back_all_tf <- ordered_back_all_tf[-1,]
      } else {
        ordered_back_all_tf <- ordered_back_all_tf$result
        ordered_back_all_tf <- ordered_back_all_tf[ordered_back_all_tf$term_size > 15 & ordered_back_all_tf$term_size < 5000 & ordered_back_all_tf$intersection_size > 2,]
      }  
      
    }

    if(saveWorks) {    
      
    grDevices::pdf(file = paste0(path,"/",output_directory,"/bp_reranked/",colnames(scMappR_vals)[i],"_reranked_path.pdf"))
    scMappR::plotBP(ordered_back_all)
    grDevices::dev.off()
    
    grDevices::pdf(file = paste0(path,"/",output_directory,"/tf_reranked/",colnames(scMappR_vals)[i],"_reranked_tf.pdf"))
    scMappR::make_TF_barplot(ordered_back_all_tf)
    grDevices::dev.off()
    
    grDevices::png(filename = paste0(path,"/",output_directory,"/bp_reranked/",colnames(scMappR_vals)[i],"_reranked_path.png"))
    scMappR::plotBP(ordered_back_all)
    grDevices::dev.off()
    
    grDevices::png(filename  = paste0(path,"/",output_directory,"/tf_reranked/",colnames(scMappR_vals)[i],"_reranked_tf.png"))
    scMappR::make_TF_barplot(ordered_back_all_tf)
    grDevices::dev.off()
    
    }
    
    l <- list(ranks = ranks, reranked_improved_genes = rownames(ranks_improve), BP_reranked = ordered_back_all,TF_reranked = ordered_back_all_tf )
    
    celltype_path[[i]] <- l
  }
  names(celltype_path) <- colnames(scMappR_vals)

  rank_increase <- celltype_path
  
  if(saveWorks) {    
    
  save(rank_increase, file = paste0(path, "/",output_directory,"/rank_increase_genes.Rdata"))
  #save(non_rank_increase, file = paste0(path, "/",output_directory,"/non_rank_increase_genes.Rdata"))
  
  }

  
  
  
  return(list(rank_increase = rank_increase, non_rank_increase = non_rank_increase))
}
