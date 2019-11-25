#' Generate scMappR Transformed Values (STVs)
#' 
#' This function takes a count matrix, signature matrix, and DEGs before generating STVs for each cell-type.
#'
#' This function completes the cell-type contextualization in scMappR -- reranking every DEG based on their fold change, likelihood the gene is in each detected cell type, average cell-type proportion, and ratio of cell-type proportion between case and control.
#' If a gene is upregulated, then it is being controlled by control/case, otherwise it is case/control.
#' This function expects that the genes within the count file, signature matrix, and DEG_list are have the same logos
#' It first takes the most "sig_matrix_size" variable genes and completes cell-type deconvolution with it using the DeconRNA-seq package before extracting cell-types with an average of greater than 0.1% of the cell-type proportion's.
#' In then subsets the signature matrix to only investigate those cell-types for the rest of the analysis. This is because cell-types that are too rare often end with fold-changes between case and control that are unreliable.
#' Then, it takes a leave-one-out approach to cell-type deconvolution such that estimated cell-type proportions are computed for every gene where they are not in the signature or bulk matrix.
#' These values are additionally plotted so you can see which genes are contributing the most highly to estimating CT proportions in your bulk sample.
#' Once these values are computed, cell type means can cell-type ratios between cases and controls for each gene are computed. You can choose to scale your odds ratios such that everything scales to the smallest non-zero value but when using absolute odds ratios it is not reccomented. 
#' Then, for every gene, scMappR is applied, the example for upreguatled genes is here: val <- cell-preferences * cell-type_proprtion * cell-type_proportion_fold-change * sign*2^abs(gene_DE$log2fc).
#' Once these values are computed, they are returned. If printplots = true, then this function also returns boplots giving the estimated cell-type proportion with each gene removed.
#'
#'
#'
#' @rdname deconvolute_and_contextualize
#' @name deconvolute_and_contextualize
#'
#' @param ordered_back_all output of the g:ProfileR function
#' @param top_bp The number of pathways you want to plot
#' 
#' @param count_file Normalized RNA-seq count matrix where rows are gene symbols and columns are individuals. Either the object tself of the path of a TSV file
#' @param signature_matrix Signature matrix (reccommended odds ratios) of cell-type specificity of genes. Either the object itself or a pathway to an RData file containing an object named "wilcoxon_rank_mat_or" -- generally internal
#' @param DEG_list An object with the first column as gene symbols within the bulk dataset (doesn't have to be in signature matrix), second column is the adjusted P-value, and the third the log2FC. Path to a tsv file containing this info is also acceptable.
#' @param case_grep. Tag in the column name for cases (i.e. samples representing upregulated) OR an index of cases.
#' @param control_grep Tag in the column name for control (i.e. samples representing downregulated) OR an index of cases.
#' @param max_proportion_change Maximum cell-type proportion change -- may be useful if there are many rare cell-types.
#' @param print_plots Whether boxplots of the estimated CT proportion for the leave-one-out method of CT deconvolution should be printed.
#' @param plot_names If plots are being printed, the prefix of their pdf files.
#' @param theSpecies -9 if using a precomputed count matrix from scMappR, human otherwise. removes ensembl symbols if appended.
#' @param make_scale convert the lowest odds ratio to 1 and scales accordingly -- strongly not reccomended and will produce warning if used
#' @param FC_coef Making scMappR Transformed Values STV's based on fold-change (TRUE) or rank := (-log10(Pval)) (FALSE) rank. After testing, we strongly reccomend to keep true.
#' @param sig_matrix_size Number of genes in signature matrix for cell-type deconvolution.
#' @param sig_distort Exponential change of odds ratios. Strongly not reccomended and will produce warnings if changed from default.
#' @param drop_unkown_celltype Whether or not to remove "unknown" cell-types from the signature matrix.
#' 
#' 
#' @return \code{deconvolute_and_contextualize} ScMappR transformed Values for every gene in every cell-type, cell-type composition with, allgenes included, average gene expression of each cell-type usng leave one out approach for each gene, and the processed signature matrix. Optional: boxplots of estimated CT proportions for each gene using a leave-one-out method \cr
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
#' load("~/scMappR/data/PBMC_scMappR_and_pathway_analysis_example.rda")
#' #data(PBMC_scMappR)
#' bulk_DE_cors <- PBMC_example$bulk_DE_cors
#' bulk_normalized <- PBMC_example$bulk_normalized
#' odds_ratio_in <- PBMC_example$odds_ratio_in
#' case_grep <- "_female"
#' control_grep <- "_male"
#' max_proportion_change <- 10
#' print_plots <- FALSE
#' theSpecies <- "human"
#'
#' norm <- deconvolute_and_contextualize(bulk_normalized, odds_ratio_in, bulk_DE_cors, case_grep = case_grep, control_grep = control_grep,
#'                                      max_proportion_change = max_proportion_change, print_plots = print_plots, theSpecies = theSpecies)
#'  }
#' @export
NULL
#' @rdname deconvolute_and_contextualize
#' @export
#' 
deconvolute_and_contextualize <- function(count_file,signature_matrix, DEG_list, case_grep, control_grep, max_proportion_change = -9, print_plots=T, plot_names="scMappR",theSpecies = "human", make_scale = F, FC_coef = T, sig_matrix_size = 3000, sig_distort = 1, drop_unkown_celltype = TRUE) {
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
  # making scMappR Transformed Values STV's based on rank (-log10(Pval)) instead of fold change
  #sig_matrix_size 
  # number of genes in signature matrix for cell-type deconvolution
  #sig_distort
  # exponential change of odds ratios. not reccomended and will produce warnings if changed from default
  # drop_unkown_celltype
  # whether or not to remove "unknown" cell-types from the signature matrix
  # Returns:
  # ScMappR transformed Values for every gene in every cell-type, cell-type composition with, allgenes included, average gene expression of each cell-type usng leave one out approach for each gene, and the processed signature matrix
  #optional: boxplots of estimated CT proportions for each gene using a leave-one-out method
  
  # load required packages

  
  rowVars <- function(x) return(apply(x, 1, var)) # get variance of rows, used later
  # load in normalized count matrices, signature matrix, and DEGs
  if(class(count_file) == "character") {
    norm_counts_i <- read.table(count_file, header = T, as.is = T, sep = "\t")
  } else {
    norm_counts_i <- count_file
  }
  if(class(signature_matrix) == "character") {
    load(signature_matrix)
    
  } else {
    wilcoxon_rank_mat_or <- signature_matrix
  }
  if(class(DEG_list) == "character") {
    DEGs <- read.table(DEG_list, header = F, as.is = T, sep = "\t")
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
    print("Rownames of signature matrix")
    print(head(RN_2))
    print("Rownames of count matrix")
    print(head(norm_counts_i))
    stop("There is no overlap between the signature and count matrix. Please make sure that gene symbols are row names and they are gene symbols of te same species.")
  }
  rownames(wilcoxon_rank_mat_or) <- RN_2
  # Removing unknown cell-types
  unknown <- grep("unknown",colnames(wilcoxon_rank_mat_or))
  if(length(unknown) > 0 & drop_unkown_celltype == TRUE) {
    print("Removing unknown cell-types")
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
  if(nrow(wilcox_or) > sig_matrix_size) {
    print(paste0("For deconvolution, we're using the top ", sig_matrix_size," most vairable signatures"))
    RVar <- matrixStats::rowVars(wilcox_var)
    wilcox_or_var <- wilcox_or[order(RVar, decreasing = T), ]
    wilcox_or_signature <- wilcox_or_var[1:sig_matrix_size,]
  } else {
    wilcox_or_signature <- wilcox_or
  }
  if(sig_distort != 1) {
    warning("Changed signature matrix distortion from 1 which is not reccomended.")
    wilcox_or_signature <- wilcox_or_signature^sig_distort
  }
  if(length(unique(colnames(wilcox_or_signature))) < length(colnames(wilcox_or_signature))) {
    warning("cell-type naming is not unique, appending unique identifier (1:ncol(signature))")
    colnames(wilcox_or_signature) <- paste0(colnames(wilcox_or_signature), "_", 1:ncol(wilcox_or_signature))
  }
  if(class(norm_counts_i) == "matrix") norm_counts_i <- as.data.frame(norm_counts_i)
  
  # cell-type deconvolution with all genes included
  all_genes_in <- DeconRNASeq::DeconRNASeq(norm_counts_i, wilcox_or_signature, fig = F)
  
  proportions <- all_genes_in$out.all
  rownames(proportions) <- colnames(norm_counts_i)
  
  propmeans <- colMeans(proportions)
  
  # keep cell-type of genes with > 0.1% of the population
  proportions <- proportions[,colMeans(proportions) > 0.001]
  print("your bulk data contains the following cell types")
  print(colnames(proportions), quote = F)
  #convert to correct datatypes for downstream analysis
  wilcox_or <- wilcox_or[,colnames(proportions)]
  wilcox_or_df <- as.data.frame(wilcox_or)
  wilcox_or_signature <- as.data.frame(wilcox_or_signature[,colnames(proportions)])
  bulk_in <- norm_counts_i
  
  DEGs <- DEGs[complete.cases(DEGs),]
  
  
  genesIn <- DEGs$gene_name[DEGs$gene_name %in% rownames(bulk_in)]
  print(length(genesIn))
  
  if(length(genesIn) == 0) {
    stop("None of your DEGs overlap with genes in your count matrix")
  }
  
  if(class(bulk_in) != "data.frame") bulk_in <- as.data.frame(bulk_in)
  if(class(wilcox_or) != "data.frame") wilcox_or <- as.data.frame(wilcox_or)
  
  
  deconvolute_gene_removed <- function(x, bulk = bulk_in, signature = wilcox_or_signature) {
    #Internal: Given an inputted gene, this function will complete DeconRNAseq witith that gene removed from the bulk and signature matrices
    # Args:
    # x: the gene symbol to be removed
    # bulk: the bulk RNA-seq dataset
    # signature: the signature matrix
    #Returns: estimated CT propotions with that gene removed
    if(x %in% rownames(bulk)) { # remove from bulk
      bulk_rem <- bulk[-which(rownames(bulk) == x),]
    } else {
      warning(paste0(x, " is not in your bulk matrix, remove gene and rerun or check gene symbol matching in general."))
      bulk_rem <- bulk
    }
    if(x %in% rownames(signature)) { # remove from signature
      signature_rem <- signature[-which(rownames(signature) == x),]
    } else {
      signature_rem <- signature
    }
    
    
    testMine <- DeconRNASeq::DeconRNASeq(bulk_rem, signature_rem, fig = F)
    
    proportions <- testMine$out.all # get proprotions
    
    print(x)
    return(proportions)
  }
  
  # run this with all genes
  iterated <- lapply(genesIn, deconvolute_gene_removed)
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
    cases_Med <- matrixStats::colMedians(tester[cases,])
    control_Med <- matrixStats::colMedians(tester[control,])
    
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
  
  print("Leave one out cell proporitons: ")
  iterated_pull <- lapply(iterated, proportion_pull)
  pull_means <- function(x) return(x$Mean)
  pull_fc <- function(x) return(x$FC) 
  # pull the means and fold-changes and bind it into a matrix
  means <- S4Vectors::do.call("rbind",lapply(iterated_pull, pull_means))
  fold_changes <- S4Vectors::do.call("rbind",lapply(iterated_pull, pull_fc))
  if( max_proportion_change != -9) { # if there is a maximum of cell-type proprotion changes, cap it at your max
    print("converting maximum CT proportion change to have a maximum odds-ratio of")
    print(max_proportion_change)
    fold_changes[fold_changes > max_proportion_change] <- max_proportion_change
    fold_changes[fold_changes < (1/max_proportion_change)] <- (1/max_proportion_change)
  }
  
  # get the average cell-type proportion for the leave one out approach with every gene removed
  cmeaned <- lapply(iterated, colMeans) 
  cmeaned_stacked <- S4Vectors::do.call("rbind", cmeaned)
  n <- colnames(cmeaned_stacked)
  print("Done")
  cmeaned_no0 <- as.data.frame(cmeaned_stacked[,colSums(cmeaned_stacked) > 0 ])
  # again, keep cell-types that have proportions
  if(length(colnames(cmeaned_stacked)) != length(unique(colnames(cmeaned_stacked)))) {
    # add a unique identifier for each cell-type if there are multiple cell-types with the same name
    print("adding code for non-unique cell-types")
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
  if(length(dup_gene_names) > 0) {
    dup_gene_names <- DEGs$gene_name[dup_gene_names]
    warning("Duplicated gene names:")
    print(dup_gene_names, quote = F)
    warning("Some gene names are duplicated -- keeping the first duplicate in the list. Check if there should be duplicate gene symbols in your dataset.")
  }
  DEGs <- DEGs[!duplicated(DEGs$gene_name),]
  rownames(DEGs) <- DEGs$gene_name
  intered <-  genesIn
  
  
  toInter_InGene<- DEGs[intered,]
  fc_co <- FC_coef
  # time to make the formula
  
  values_with_preferences <- function(gene, FC_coef = fc_co) {
    # Internal: Makes scMappR Transformed Variables (STVs): combines fold change of DEG, likelihood DEG is in cell-type, cell-type proportion with DEG left out, and odds-ratio of cell-type proportions 
    # Args:
    # gene: the gene symbol to generate a STV for
    # FC_coef: if the STV is based on Fold-Change (reccomended) or Rank
    # Returns:
    # STVs for every cell-type within a single gene.
    print(gene)
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
    if(up == T) { # if it's upregulated 
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
  print("Adjusting Coefficients:")
  vals_out <- lapply(toInter_InGene$gene_name, values_with_preferences)
  print("Done")
  vals_out_mat <- S4Vectors::do.call("rbind", vals_out)
  rownames(vals_out_mat) <- toInter_InGene$gene_name
  
  colnames(vals_out_mat) <- colnames(wilcox_or)
  all_reordered <-  vals_out_mat
  comp <- plot_names
  
  if(print_plots == T) {
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
        bottom3 <- tail(sort(cm_stack), 3)
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

      pdf(paste0("deconvolute_generemove_quantseq_",names,".pdf"))
      g <- ggplot2::ggplot(all_stack, ggplot2::aes(factor(cell_type), proportion)) + ggplot2::geom_boxplot()  + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8))
      # generate barplot for every cell-type combined
      plot(g)
      dev.off()
      
      for(i in unique(all_stack$cell_type)) {
        cm_one <- all_stack[all_stack$cell_type == i,]
        pdf(paste0("deconvolute_generemov_quantseq_", i,"_",names,".pdf"))
        g <- ggplot2::ggplot(cm_one, ggplot2::aes(factor(cell_type), proportion)) + ggplot2::geom_boxplot() + ggplot2::geom_text(ggplot2::aes(label=label),hjust=-0.2) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 1, size = 12))
        # generate barplot for one cell-type at a time
        plot(g)
        dev.off()
        print(i)
      }
      return("Done!")
    }
    
    boxplot_values(cmeaned_stacked, paste0(plot_names,"_average"))
    #boxplot_values(means, paste0(plot_names,"_median"))
    boxplot_values(fold_changes, paste0(plot_names,"_foldchange"))
  }
  l <- list(scMappR_transformed_values = all_reordered, 
            cellType_Proportions = proportions,
            leave_one_out_proportions = cmeaned_stacked,
            processed_signature_matrix = scaled_odds_ratio)   
  return(l)
  
  
}

