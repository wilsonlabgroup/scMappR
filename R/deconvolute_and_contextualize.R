#' Generate cell weighted Fold-Changes (cwFold-changes)
#' 
#' This function takes a count matrix, signature matrix, and differentially expressed genes (DEGs) before generating cwFold-chagnes for each cell-type.
#'
#' This function completes the cell-type contextualization in scMappR -- reranking every DEG based on their fold change, likelihood the gene is in each detected cell-type, average cell-type proportion, and ratio of cell-type proportion between case and control.
#' If a gene is upregulated then it is being controlled by control/case, otherwise it is case/control.
#' cwFold-change's are generated for genes that are in both the count matrix and in the list of DEGs. It does not have to also be in the signature matrix.
#' First, this function will estimate cell-type proportions with all genes included before estimating changes in cell-type proportion between case/control using a t-test.
#' Then, it takes a leave-one-out approach to cell-type deconvolution such that estimated cell-type proportions are computed for every inputted DEG.
#' Optionally, the differences between cell-type proprtion before and after a gene is removed is plotted in boxplots.
#' Then, for every gene, cwFold-change's are computed with the following formula (the example for upreguatled genes)
#' val <- cell-preferences * cell-type_proprtion * cell-type_proportion_fold-change * sign*2^abs(gene_DE$log2fc).
#' A matrix of cwFold-change's for all DEGs are returned.
#'
#'
#'
#' @rdname deconvolute_and_contextualize
#' @name deconvolute_and_contextualize
#'
#' @param count_file Normalized deconvolute_and_contextualize. RNA-seq count matrix where rows are gene symbols and columns are individuals. Either the object tself of the path of a .tsv file.
#' @param signature_matrix Signature matrix (recommended odds ratios) of cell-type specificity of genes. Either the object itself or a pathway to a .RData file containing an object named "wilcoxon_rank_mat_or" - generally internal.
#' @param DEG_list An object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted P-value, and the third the log2FC. Path to a tsv file containing this info is also acceptable.
#' @param case_grep Tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases.
#' @param control_grep Tag in the column name for control (i.e. samples representing downregulated) OR an index of cases.
#' @param max_proportion_change Maximum cell-type proportion change. May be useful if there are many rare cell-types.
#' @param print_plots Whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed (T/F).
#' @param plot_names If plots are being printed, the pre-fix of their .pdf files.
#' @param theSpecies -9 if using a precomputed count matrix from scMappR, human otherwise. Removes ensembl symbols if appended.
#' @param make_scale Convert the lowest odds ratio to 1 and scales accordingly -- strongly not recommended and will produce warning if used.
#' @param FC_coef Making cwFold-changes based on fold-change (TRUE) or rank := (-log10(Pval)) (FALSE) rank. After testing, we strongly recommend to keep true (T/F).
#' @param sig_matrix_size Number of genes in signature matrix for cell-type deconvolution.
#' @param sig_distort Exponential change of odds ratios. Strongly not recomended and will produce warnings if changed from default.
#' @param drop_unknown_celltype Whether or not to remove "unknown" cell-types from the signature matrix (T/F).
#' @param toSave Allow scMappR to write files in the current directory (T/F).
#' @param path If toSave == TRUE, path to the directory where files will be saved.
#' @param deconMethod Which RNA-seq deconvolution method to use to estimate cell-type proporitons. Options are "WGCNA", "DCQ", or "DeconRNAseq"
#' 
#' @return List with the following elements:
#' \item{cellWeighted_Foldchange}{data frame of cellweightedFold changes for each gene.}
#' \item{cellType_Proportions}{data frame of cell-type proportions from DeconRNA-seq.}
#' \item{leave_one_out_proportions}{data frame of average cell-type proportions for case and control when gene is removed.}
#' \item{processed_signature_matrix}{signature matrix used in final analysis.}
#'  
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_text theme coord_flip labs element_text
#' @importFrom pheatmap pheatmap
#' @importFrom graphics barplot plot
#' @importFrom Seurat AverageExpression CreateSeuratObject PercentageFeatureSet SCTransform SelectIntegrationFeatures PrepSCTIntegration FindIntegrationAnchors IntegrateData DefaultAssay RunPCA RunUMAP FindNeighbors FindClusters ScaleData FindMarkers
#' @importFrom GSVA gsva
#' @importFrom stats fisher.test median p.adjust reorder t.test sd var complete.cases ks.test
#' @importFrom utils combn read.table write.table head tail
#' @importFrom downloader download
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom gprofiler2 gost
#' @importFrom gProfileR gprofiler
#' @importFrom pcaMethods prep pca R2cum
#' @importFrom limSolve lsei
#' @importFrom pbapply pblapply
#' @importFrom ADAPTS estCellPercent
#' @importFrom reshape melt
#' 
#' @examples 
#' \donttest{
#' data(PBMC_example)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#' norm <- deconvolute_and_contextualize(bulk_normalized, odds_ratio_in, bulk_DE_cors,
#'                                     case_grep = case_grep, control_grep = control_grep,
#'                                      max_proportion_change = max_proportion_change,
#'                                       print_plots = print_plots, 
#'                                      theSpecies = theSpecies, toSave = FALSE)
#'
#' }                                      
#' @export
#' 
deconvolute_and_contextualize <- function(count_file,signature_matrix, DEG_list, case_grep, control_grep, max_proportion_change = -9, print_plots=T, plot_names="scMappR",theSpecies = "human", make_scale = FALSE, FC_coef = T, sig_matrix_size = 3000, sig_distort = 1, drop_unknown_celltype = TRUE, toSave = FALSE, path = NULL, deconMethod = "DeconRNASeq") {
  # This function completes the cell-type contextualization in scMappR -- reranking every DEG based on their fold change, likelihood the gene is in each detected cell type, average cell-type proportion, and ratio of cell-type proportion between case and control.
  # such that if a gene is upregulated, then it is being controlled by control/case, otherwise it is case/control
  # This function expects that the genes within the count file, signature matrix, and DEG_list are have the same logos
  # it first takes the most "sig_matrix_size" variable genes and completes cell-type deconvolution with it using the DeconRNA-seq package before extracting cell-types with an average of greater than 0.1% of the cell-type proportion's
  # in then subsets the signature matrix to only investigate those cell-types for the rest of the analysis. This is because cell-types that are too rare often end with fold-changes between case and control that are unreliable.
  # Then, it takes a leave-one-out approach to cell-type deconvolution such that estimated cell-type proportions are computed for every gene where they are not in the signature or bulk matrix.
  # these values are additionally plotted so you can see which genes are contributing the most highly to estimating CT proportions in your bulk sample
  # once these values are computed, cell type means can cell-type ratios between cases and controls for each gene are computed. You can choose to scale your odds ratios such that everything scales to the smallest non-zero value but when using absolute odds ratios it is not reccomented 
  # then, for every gene, scMappR is applied, the example for upreguatled genes is here: val <- cell-preferences * cell-type_proprtion * cell-type_proportion_fold-change * sign*2^abs(gene_DE$log2fc)
  # once these values are computed, they are returned
  # if printplots = true, then this function also returns boplots giving the estimated cell-type proportion with each gene removed
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
  # tag in the column name for control (i.e. samples representing downregulated) OR an index of cases
  #max_proportion_change
  # maximum cell-type proportion change -- may be useful if there are many rare cell-types
  #print_plots 
  # whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed
  #plot_names
  # if plots are being printed, the prefix of their pdf files
  #theSpecies
  # -9 if using a precomputed count matrix from scMappR, human otherwise.
  # removes ensembl symbols
  #make_scale
  # convert the lowest odds ratio to 1 and scales accordingly -- strongly not suggested and will produce warning if used
  #FC_coef
  # making cell weighted fold-changes cwFold-change's based on rank (-log10(Pval)) instead of fold change
  #sig_matrix_size 
  # number of genes in signature matrix for cell-type deconvolution
  #sig_distort
  # exponential change of odds ratios. not reccomended and will produce warnings if changed from default
  # drop_unknown_celltype
  # whether or not to remove "unknown" cell-types from the signature matrix
  # Returns:
  # cell weighted fold-changes for every gene in every cell-type, cell-type composition with, allgenes included, average gene expression of each cell-type usng leave one out approach for each gene, and the processed signature matrix
  #optional: boxplots of estimated CT proportions for each gene using a leave-one-out method
  
  # load required packages
  
  count_class <- class(count_file)[1] %in% c("character", "data.frame", "matrix")
  if((count_class[1] == FALSE )[1]) {
    stop("count_file must be of class character, data.frame, or matrix.")
  }
  signature_class <- class(signature_matrix)[1] %in% c("character", "data.frame", "matrix")
  if((signature_class[1] == FALSE)[1]) {
    stop("count_file must be of class character, data.frame, or matrix.")
  }
  DEG_list_class <- class(DEG_list)[1] %in% c("character", "data.frame", "matrix") 
  if((DEG_list_class[1] == FALSE)[1]) {
    stop("DEG_list must be of class character, data.frame, or matrix.")
  }
  case_grep_class <- class(case_grep)[1] %in% c("character", "numeric", "integer")
  if((case_grep_class[1] == FALSE)[1]) {
    stop("case_grep must be of class character (as a single character designating cases in column names) or of class numeric (integer matrix giving indeces of cases).")
  }
  control_grep_class <- class(control_grep)[1] %in% c("character", "numeric", "integer")
  
  if((control_grep_class == FALSE)[1]) {
    stop("control_grep must be of class character (as a single character designating controls in column names) or of class numeric (integer matrix giving indeces of controls).")
  }
  
  
  if((all(is.numeric(sig_distort), is.numeric(max_proportion_change), is.numeric(sig_matrix_size))[1] == FALSE)[1]) {
    stop("sig_distort, max_proportion_change, and sig_matrix_size must all be of class numeric" )
  }
  if((all(is.logical(print_plots), is.logical(make_scale), is.logical(FC_coef), is.logical(drop_unknown_celltype), is.logical(toSave))[1] == FALSE)[1]) {
    stop("print_plots, make_scale, FC_coef, drop_unknown_celltype, toSave must all be of class logical." )
  }
  
  if(is.character(plot_names)[1] == FALSE) {
    stop("plot_names must be of class character.")
  }
  if((!(theSpecies %in% c("human", "mouse")))[1]) {
    if((theSpecies != -9)) {
      stop("species_name is not 'human' 'mouse' or '-9' (case sensitive), please try again with this filled.")
    }
  }
  
  if(toSave == TRUE) {
    if(is.null(path)) {
      stop("scMappR is given write permission by setting toSave = TRUE but no directory has been selected (path = NULL). Pick a directory or set path to './' for current working directory")
    }
    if(!dir.exists(path)) {
      stop("The selected directory does not seem to exist, please check set path.")
    }
  }
  
  rowVars <- function(x) return(apply(x, 1, stats::var)) # get variance of rows, used later
  colMedians <- function(x) return(apply(x, 2, stats::median)) # medians of cols, used later
  colVars <- function(x) return(apply(x, 2 , stats::var))
  # load in normalized count matrices, signature matrix, and 0
  
  if(is.character(count_file)) {
    warning("reading count file table, assuming first column is the gene symbols -- not reccomended.")
    norm_counts_i <- utils::read.table(count_file, header = T, as.is = T, sep = "\t")
    rownames(norm_counts_i) <- norm_counts_i[,1]
    norm_counts_i <- norm_counts_i[,-1]
  } else {
    norm_counts_i <- count_file
  }
  if(is.character(signature_matrix)) {
    warning("loading in signature matrix from Rdata file -- not reccomeneded.")
    load(signature_matrix)
    
  } else {
    wilcoxon_rank_mat_or <- signature_matrix
  }
  if(is.character(DEG_list)) {
    DEGs <- utils::read.table(DEG_list, header = FALSE, as.is = T, sep = "\t")
  } else {
    DEGs <- as.data.frame(DEG_list)
  }
  
  colnames(DEGs) <- c("gene_name", "padj", "log2fc")
  sm <- wilcoxon_rank_mat_or
  
  # If internal, get gene name and species
  if(theSpecies == -9) {
    gsym <- get_gene_symbol(sm)
    RN_2 <- gsym$rowname
    theSpecies <- gsym$species
  } else {
    RN_2 <- rownames(sm)
  }
  toInter <- intersect(RN_2, rownames(norm_counts_i))
  
  if(length(toInter)==0) { 
    # if there isn't overlap in the genes in the signature matrices and count matrix
    message("Rownames of signature matrix")
    message(utils::head(RN_2))
    message("Rownames of count matrix")
    message(utils::head(norm_counts_i))
    stop("There is no overlap between the signature and count matrix. Please make sure that gene symbols are row names and they are gene symbols of te same species.")
  }
  rownames(wilcoxon_rank_mat_or) <- RN_2
  # Removing unknown cell-types
  unknown <- grep("unknown",colnames(wilcoxon_rank_mat_or))
  if((length(unknown) > 0 & drop_unknown_celltype == TRUE)[1]) {
    message("Removing unknown cell-types")
    wilcox_or <- wilcoxon_rank_mat_or[,-unknown]
  } else {
    wilcox_or <- wilcoxon_rank_mat_or
  }
  wilcox_var <- wilcox_or
  wilcox_or <- as.data.frame(wilcox_or)
  ############
  #### Make two wilcox_or matrices -- one with the top "3,000" most variable genes and the other full one 
  #############
  wilcox_or[wilcox_or < 0 ] <- 0
  if((nrow(wilcox_or) > sig_matrix_size)[1]) {
    message(paste0("For deconvolution, we're using the top ", sig_matrix_size," most vairable signatures"))
    RVar <- rowVars(wilcox_var)
    wilcox_or_var <- wilcox_or[order(RVar, decreasing = T), ]
    wilcox_or_signature <- wilcox_or_var[1:sig_matrix_size,]
  } else {
    wilcox_or_signature <- wilcox_or
  }
  if(sig_distort != 1) {
    warning("Changed signature matrix distortion from 1 which is not reccomended.")
    wilcox_or_signature <- wilcox_or_signature^sig_distort
  }
  if((length(unique(colnames(wilcox_or_signature))) < length(colnames(wilcox_or_signature)))[1]) {
    warning("cell-type naming is not unique, appending unique identifier (1:ncol(signature))")
    colnames(wilcox_or_signature) <- paste0(colnames(wilcox_or_signature), "_", 1:ncol(wilcox_or_signature))
  }
  if(is.matrix(norm_counts_i)) norm_counts_i <- as.data.frame(norm_counts_i)
  
  cVar_wilcox <- colVars(wilcox_or_signature)
  wilcox_or_signature <-   wilcox_or_signature[,cVar_wilcox > 0]
  # cell-type deconvolution with all genes included
  
  if(deconMethod == "DCQ") {
    all_genes_in <- t(ADAPTS::estCellPercent(wilcox_or_signature,norm_counts_i, method = "DCQ") / 100)
    all_genes_in <- all_genes_in[,colnames(all_genes_in) != "others"]
    
  }
  if(deconMethod == "DeconRNASeq") {
    all_genes_in <- DeconRNAseq_CRAN(as.data.frame(norm_counts_i),as.data.frame(wilcox_or_signature))$out.all
    
  }
  if(deconMethod == "WGCNA") {
    all_genes_in <- t(ADAPTS::estCellPercent(wilcox_or_signature,norm_counts_i, method = "proportionsInAdmixture") / 100)
    all_genes_in <- all_genes_in[,colnames(all_genes_in) != "others"]
    
  }
  #Decon <- ADAPTS::estCellPercent(bulk_normalized,odds_ratio_in, method = "DeconRNASeq")
  
  #all_genes_in <- DeconRNAseq_CRAN(norm_counts_i, wilcox_or_signature)
  
  proportions <- all_genes_in
  
  rownames(proportions) <- colnames(norm_counts_i)
  
  propmeans <- colMeans(proportions)
  
  # keep cell-type of genes with > 0.1% of the population
  proportions <- proportions[,colMeans(proportions) > 0.001]
  
  message("your bulk data contains the following cell types")
  message(paste(colnames(proportions), " "))
  #convert to correct datatypes for downstream analysis
  wilcox_or <- wilcox_or[,colnames(proportions)]
  wilcox_or_df <- as.data.frame(wilcox_or)
  wilcox_or_signature <- as.data.frame(wilcox_or_signature[,colnames(proportions)])
  bulk_in <- norm_counts_i
  
  DEGs <- DEGs[stats::complete.cases(DEGs),]
  
  
  genesIn <- DEGs$gene_name[DEGs$gene_name %in% rownames(bulk_in)]
  message(length(genesIn))
  
  if(length(genesIn) == 0) {
    stop("None of your DEGs overlap with genes in your count matrix")
  }
  
  if(!is.data.frame(bulk_in)) bulk_in <- as.data.frame(bulk_in)
  if(!is.data.frame(wilcox_or)) wilcox_or <- as.data.frame(wilcox_or)
  max_bulk <- max(bulk_in)
  if(max_bulk < 50) {
    warning("You most highly expressed gene is < 50 suggesting your counts are log2 normalized. Removing log2 transformation")
    message("You most highly expressed gene is < 50 suggesting your counts are log2 normalized. Removing log2 transformation. We suggest you input non-log transformed count data.")
    bulk_in <- 2^bulk_in
    
  }
  
  deconvolute_gene_removed <- function(x, bulk = bulk_in, signature = wilcox_or_signature) {
    #Internal: Given an inputted gene, this function will complete DeconRNASeq witith that gene removed from the bulk and signature matrices
    # Args:
    # x: the gene symbol to be removed
    # bulk: the bulk RNA-seq dataset
    # signature: the signature matrix
    #Returns: estimated CT propotions with that gene removed
    #message(x)
    if(x %in% rownames(bulk)) { # remove from bulk
      
      bulk_rem <- bulk[-which(rownames(bulk) == x),]
    } else {
      warning(paste0(x, " is not in your bulk matrix, remove gene and rerun or check gene symbol matching in general."))
      bulk_rem <- bulk
    }
    if(x %in% rownames(signature)) { # remove from signature
      signature_rem <- signature[-which(rownames(signature) == x),]
      signature_var <- apply(signature_rem, 2, stats::var)
      if(any(signature_var == 0)[1]) { 
        # if the DEG we removed is the only cell-type marker
        signature_rem <- signature
        bulk_rem <- bulk
      }
    } else {
      signature_rem <- signature
    }
    
    ##
    if(deconMethod == "DCQ") {
      proportions <- t(ADAPTS::estCellPercent(signature_rem,bulk_rem, method = "DCQ") / 100)
      proportions <- proportions[,colnames(proportions) != "others"]
      
    }
    if(deconMethod == "DeconRNASeq") {
      proportions <- DeconRNAseq_CRAN(bulk_rem,signature_rem)$out.all
      
    }
    if(deconMethod == "WGCNA") {
      proportions <- t(ADAPTS::estCellPercent(signature_rem,bulk_rem, method = "proportionsInAdmixture") / 100)
      proportions <- proportions[,colnames(proportions) != "others"]
      
    }
    ##
    #testMine <- DeconRNAseq_CRAN(bulk_rem, signature_rem)
    
    #proportions <- testMine$out.all # get proprotions
    
    #message(x)
    return(proportions)
  }
  
  # run this with all genes
  iterated <- pbapply::pblapply(genesIn, deconvolute_gene_removed)
  names(iterated) <- genesIn
  proportion_pull <- function(tester) {
    #Internal: given a list of cell-type proprotions using a leave-one-out approach, pull the average
    # cell-type proportions in case/control, average overall, and the odds ratio that the gene is found in cases
    # Args: 
    # tester: a lis of cell-type proporitons for each leave one out for each gene
    # Returns: 
    # cell-type proportions for ever gene given case/control and their odds ratio
    
    rownames(tester) <- colnames(bulk_in)
    cases <- grep(case_grep, rownames(tester))
    control <- grep(control_grep, rownames(tester))
    if((any(length(cases) < 2, length(control) < 2))[1]) {
      stop("There is fewer than two cases or controls, please check 'case_grep' or 'control_grep'.")
    }
    cases_Med <- colMedians(tester[cases,])
    control_Med <- colMedians(tester[control,])
    
    # If the cell-type is not detected in case/control at all, replace with the lowest value overall
    # this removes the chance of an infinite fold-change in CT proportion
    M <- min(c(min(cases_Med[cases_Med != 0]), min(control_Med[control_Med != 0])))
    
    cases_Med[cases_Med == 0] <- M
    control_Med[control_Med == 0] <- M
    
    FC <- cases_Med/control_Med
    Mean_CC <- (cases_Med + control_Med) / 2
    l <- list(cases = cases_Med, control = control_Med, FC = FC, Mean = Mean_CC)
    return(l) 
  }
  
  message("Leave one out cell proporitons: ")
  iterated_pull <- pbapply::pblapply(iterated, proportion_pull)
  pull_means <- function(x) return(x$Mean)
  pull_fc <- function(x) return(x$FC) 
  # pull the means and fold-changes and bind it into a matrix
  means <- do.call("rbind",lapply(iterated_pull, pull_means))
  fold_changes <- do.call("rbind",lapply(iterated_pull, pull_fc))
  if( max_proportion_change != -9) { # if there is a maximum of cell-type proprotion changes, cap it at your max
    message("converting maximum CT proportion change to have a maximum odds-ratio of")
    message(max_proportion_change)
    fold_changes[fold_changes > max_proportion_change] <- max_proportion_change
    fold_changes[fold_changes < (1/max_proportion_change)] <- (1/max_proportion_change)
  }
  
  # get the average cell-type proportion for the leave one out approach with every gene removed
  cmeaned <- lapply(iterated, colMeans) 
  cmeaned_stacked <- do.call("rbind", cmeaned)
  n <- colnames(cmeaned_stacked)
  message("Done")
  cmeaned_no0 <- as.data.frame(cmeaned_stacked[,colSums(cmeaned_stacked) > 0 ])
  # again, keep cell-types that have proportions
  if((length(colnames(cmeaned_stacked)) != length(unique(colnames(cmeaned_stacked))))[1]) {
    # add a unique identifier for each cell-type if there are multiple cell-types with the same name
    message("adding code for non-unique cell-types")
    colnames(fold_changes) <- colnames(wilcox_or) <- colnames(means) <- colnames(cmeaned_stacked) <- paste0(n,"_",1:ncol(cmeaned_stacked))
  } else {
    colnames(fold_changes) <- colnames(wilcox_or) <- colnames(means) <- colnames(cmeaned_stacked) <- n
  }
  if(make_scale == TRUE) { 
    warning("Scaling all odds ratios to the minimum odds ratio -- not reccomended")
    scaled <- min(wilcox_or[wilcox_or != 0])
    
    scaled_odds_ratio <- wilcox_or/scaled
  } else {
    scaled_odds_ratio <- wilcox_or
  }
  # remove duplicated DEGs
  dup_gene_names <- which(duplicated(DEGs$gene_name))
  if((length(dup_gene_names) > 0)[1]) {
    dup_gene_names <- DEGs$gene_name[dup_gene_names]
    warning("Duplicated gene names:")
    message(paste(dup_gene_names, " "))
    warning("Some gene names are duplicated -- keeping the first duplicate in the list. Check if there should be duplicate gene symbols in your dataset.")
  }
  DEGs <- DEGs[!duplicated(DEGs$gene_name),]
  rownames(DEGs) <- DEGs$gene_name
  intered <-  genesIn
  
  
  toInter_InGene<- DEGs[intered,]
  fc_co <- FC_coef
  # time to make the formula
  
  values_with_preferences <- function(gene, FC_coef = fc_co) {
    # Internal: Makes scMappR Transformed Variables (cwFold-changes): combines fold change of DEG, likelihood DEG is in cell-type, cell-type proportion with DEG left out, and odds-ratio of cell-type proportions 
    # Args:
    # gene: the gene symbol to generate a cwFold-change for
    # FC_coef: if the cwFold-change is based on Fold-Change (reccomended) or Rank
    # Returns:
    # cwFold-changes for every cell-type within a single gene.
    #message(gene)
    if(gene %in% rownames(scaled_odds_ratio)) {
      scaled_pref <- scaled_odds_ratio[gene,] # extract gene from signature
    } else {
      scaled_pref <- rep(1, ncol(scaled_odds_ratio)) # if it's not in signautre assume that it's equally likely to be found in each cell-type
    }
    gene_DE <- DEGs[gene,] # DE value of gene
    prop <- means[gene,] # cell-type proportion average of gene
    prop_fc <- fold_changes[gene, ] # relative ratio of case/control in gene
    up <- gene_DE$log2fc > 0 # if gene is upregualted 
    sign <- -1
    if((up == T)[1]) { # if it's upregulated 
      prop_fc <- 1/prop_fc # flip so the ratio is control/case
      sign <- 1
    } 
    if(FC_coef == TRUE) {
      
      val <- scaled_pref * prop * prop_fc * sign*2^abs(gene_DE$log2fc)
    } else {
      warning("using P-FDR instead of fold-change -- Not reccomended")
      val <- scaled_pref * prop * prop_fc * sign*-1*log10(gene_DE$padj)  
    }
    return(val)
  }
  message("Adjusting Coefficients:")
  vals_out <- pbapply::pblapply(toInter_InGene$gene_name, values_with_preferences)
  message("Done")
  vals_out_mat <- do.call("rbind", vals_out)
  rownames(vals_out_mat) <- toInter_InGene$gene_name
  
  colnames(vals_out_mat) <- colnames(wilcox_or)
  all_reordered <-  vals_out_mat
  comp <- plot_names
  printing <- print_plots == TRUE & toSave == TRUE
  if(printing[1] == TRUE) {
    # If you want to print arplots
    boxplot_values <- function(cmeaned_stacked, names) {
      #Internal: get the values of cell-type proportions with a leave-one-out method and trint in boxplot
      # print name of genes with top and bottom 3 DEGs
      # Args: 
      # cmeaned_stacked: average cell-type propotions across all samples using leave-one-out approach
      # names: the prefix of the pdf 
      # Returns: 
      # Boxplots for estimated cell-type proportions for a leave-one-out cell-type
      # for every cell-type put together and for one cell-type individually
      all_stack <- c()
      for(i in 1:ncol(cmeaned_stacked)) {
        cm_stack <- cmeaned_stacked[,i]
        # name top 3 genes, replace non top-3 with ""
        top3 <- head(sort(cm_stack), 3)
        bottom3 <- utils::tail(sort(cm_stack), 3)
        empty <- rep("", length(cm_stack))
        names(empty) <- names(cm_stack)
        empty[names(top3)] <- names(top3)
        empty[names(bottom3)] <- names(bottom3)
        y <- cbind(names(cm_stack), unname(cm_stack), colnames(cmeaned_stacked)[i], unname(empty))
        all_stack <- rbind(all_stack, y)
        # combine into ggplot2 compatible barchart
      }
      all_stack <- as.data.frame(all_stack)
      # make sure that data in the barplots are all of the right datatypr
      colnames(all_stack) <- c("gene", "proportion", "cell_type", "label")
      all_stack$gene <- tochr(all_stack$gene)
      all_stack$proportion <- toNum(all_stack$proportion)
      all_stack$cell_type <- tochr(all_stack$cell_type)
      all_stack$label <- tochr(all_stack$label)
      cell_type <- all_stack$cell_type
      proportion<- all_stack$proportion
      label <- all_stack$label
      grDevices::pdf(paste0(path,"/","deconvolute_generemove_quantseq_",names,".pdf"))
      g <- ggplot2::ggplot(all_stack, ggplot2::aes(factor(cell_type), proportion)) + ggplot2::geom_boxplot()  + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8))
      # generate barplot for every cell-type combined
      print(g)
      grDevices::dev.off()
      
      for(i in unique(all_stack$cell_type)) {
        cm_one <- all_stack[all_stack$cell_type == i,]
        cell_type <- cm_one$cell_type
        proportion <- cm_one$proportion
        label <- cm_one$label
        g <- ggplot2::ggplot(cm_one, ggplot2::aes(factor(cell_type), proportion)) + ggplot2::geom_boxplot() + ggplot2::geom_text(ggplot2::aes(label=label),hjust=-0.2) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 1, size = 12))
        
        grDevices::pdf(paste0(path,"/","deconvolute_generemov_quantseq_", i,"_",names,".pdf"))
        # generate barplot for one cell-type at a time
        print(g)
        grDevices::dev.off()
        message(i)
      }
      return("Done!")
    }
    
  }
  l <- list(cellWeighted_Foldchange = all_reordered, 
            cellType_Proportions = proportions,
            leave_one_out_proportions = cmeaned_stacked,
            processed_signature_matrix = scaled_odds_ratio)   
  return(l)
  
  
}


