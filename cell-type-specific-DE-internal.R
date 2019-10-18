A <- Sys.time()
tochr <- function(x) {
  if(class(x) == "character") return(x)
  if(class(x) == "factor") return(as.character(levels(x))[x])
  if(class(x) == "numeric") return(as.character(x))
} 

toNum <- function(x) {
  if(class(x) == "character") return(as.numeric(x))
  if(class(x) == "factor") return(as.numeric(levels(x))[x])
  if(class(x) == "numeric") return(x)
} 

pull_means <- function(x) return(x$Mean)
pull_fc <- function(x) return(x$FC) # Temporary bug -- There's an issue with fold-changes where rare cell-types get huge fold-changes from case/control


deconvolute_and_contextualize <- function(count_file,odds_ratio, DEG_list, case_grep, control_grep, max_proportion_change = -9, print_plots=T, plot_names="scMappR",theSpecies = -9) {
  if(class(count_file) == "character") {
    norm_counts_i <- read.table(count_file, header = T, as.is = T, sep = "\t")
  } else {
    norm_counts_i <- count_file
  }
  if(class(odds_ratio) == "character") {
    load(odds_ratio)
    
  } else {
    wilcoxon_map_rnk_or <- odds_ratio
  }
  if(class(DEG_list) == "character") {
    DEGs <- read.table(DEG_list, header = F, as.is = T, sep = "\t")
  } else {
    DEGs <- as.data.frame(DEG_list)
  }
  
  colnames(DEGs) <- c("gene_name", "padj", "log2fc")
  sm <- wilcoxon_map_rnk_or
  if(theSpecies == -9) {
    # if it was a panglao dataset, remove ensembl gene name
    # and based on the type of the ensembl name, infer species
    the_human <- length(grep("ENSG00", rownames(sm)))
    the_mouse <- length(grep("ENSMUSG00", rownames(sm)))
    if(the_human >= the_mouse) {
      theSpecies <- "human"
    }
    if(the_mouse > the_human) {
      theSpecies <- "mouse"
    }
    # make sure that mitochondrial genes are flagged
    RN = rownames(sm)
    RN_1 <- sub('(.*)[.](.*)','\\1',RN)
    RN_2 <- substr(RN_1, 1, nchar(RN_1) - 19)
    rownames(sm) <- RN_2
    
  }
  RN_2 <- rownames(sm)
  rownames(wilcoxon_map_rnk_or) <- RN_2
  wilcox_or <- wilcoxon_map_rnk_or[,which(colnames(wilcoxon_map_rnk_or) != "unknown")]
  wilcox_or <- as.data.frame(wilcox_or)
  wilcox_or[wilcox_or < 0 ] <- 0
  if(class(norm_counts_i) == "matrix") norm_counts_i <- as.data.frame(norm_counts_i)
  all_genes_in <- DeconRNASeq(norm_counts_i, wilcox_or)
  
  proportions <- all_genes_in$out.all
  rownames(proportions) <- colnames(norm_counts_i)
  
  propmeans <- colMeans(proportions)
  proportions <- proportions[,colSums(proportions) > 0]
  print("your bulk data contains the following cell types")
  print(colnames(proportions), quote = F)
  wilcox_or <- wilcox_or[,colnames(proportions)]
  wilcox_or_df <- as.data.frame(wilcox_or)
  bulk_in <- norm_counts_i
  print(head(DEGs$gene_name))
  print(head(rownames(wilcox_or)))
  DEGs <- DEGs[complete.cases(DEGs),]
  genesIn <- intersect(DEGs$gene_name, rownames(wilcox_or))
  print(length(genesIn))
  
  if(length(genesIn) < ncol(wilcox_or)) {
    print("there are fewer DEGs")
    print(length(genesIn))
    print("than there are cell-types")
    print(ncol(wilcox_or))
    print("please remove cell types, lower DEG threshold, or select 'all genes' option")
    print("dustin -- build the all genes option before you forget... dumbass")
    #return("more DEGS pls")
  }
  genesOut <- DEGs$gene_name[!(DEGs$gene_name %in% rownames(wilcox_or))]
  if(class(bulk_in) != "data.frame") bulk_in <- as.data.frame(bulk_in)
  if(class(wilcox_or) != "data.frame") wilcox_or <- as.data.frame(wilcox_or)
  deconvolute_gene_removed <- function(x, bulk = bulk_in, signature = wilcox_or) {
    print(x)
    bulk_rem <- bulk[-which(rownames(bulk) == x),]
    if(nrow(bulk_rem) == 0) {
      print("Warning: the following gene is not in your bulk count matrix")
      bulk_rem <- bulk
    }
    signature_rem <- signature[-which(rownames(signature) == x),]
    testMine <- DeconRNASeq(bulk_rem, signature_rem)
    proportions <- testMine$out.all
    print(x)
    return(proportions)
  }
  print("a")
  iterated <- lapply(genesIn, deconvolute_gene_removed)
  print("B")
  names(iterated) <- genesIn
  proportion_pull <- function(tester) {
    rownames(tester) <- colnames(bulk_in)
    cases <- grep(case_grep, rownames(tester))
    control <- grep(control_grep, rownames(tester))
    cases_Med <- colMedians(tester[cases,])
    control_Med <- colMedians(tester[control,])
    
    M <- min(c(min(cases_Med[cases_Med != 0]), min(control_Med[control_Med != 0])))
    
    cases_Med[cases_Med == 0] <- M
    control_Med[control_Med == 0] <- M
    
    FC <- cases_Med/control_Med
    Mean_CC <- (cases_Med + control_Med) / 2
    l <- list(cases = cases_Med, control = control_Med, FC = FC, Mean = Mean_CC)
    return(l) 
  }
  
  
  iterated_pull <- lapply(iterated, proportion_pull)
  pull_means <- function(x) return(x$Mean)
  pull_fc <- function(x) return(x$FC) 
  print(iterated_pull)
  means <- do.call("rbind",lapply(iterated_pull, pull_means))
  fold_changes <- do.call("rbind",lapply(iterated_pull, pull_fc))
  print(fold_changes)
  if( max_proportion_change != -9) {
    fold_changes[fold_changes > max_proportion_change] <- max_proportion_change
  }
  cmeaned <- lapply(iterated, colMeans)
  cmeaned_stacked <- do.call("rbind", cmeaned)
  n <- colnames(cmeaned_stacked)
  print("c")
  cmeaned_no0 <- as.data.frame(cmeaned_stacked[,colSums(cmeaned_stacked) > 0 ])
  
  colnames(fold_changes) <- colnames(wilcox_or) <- colnames(means) <- colnames(cmeaned_stacked) <- paste0(n,"_",1:ncol(cmeaned_stacked))
  
  scaled <- min(wilcox_or[wilcox_or != 0])
  
  scaled_odds_ratio <- wilcox_or/scaled
  
  DEGs <- DEGs[!duplicated(DEGs$gene_name),]
  rownames(DEGs) <- DEGs$gene_name
  
  intered <- intersect(rownames(scaled_odds_ratio), DEGs$gene_name )
  
  toInter_OR <- scaled_odds_ratio[intered,]
  toInter_InGene<- DEGs[intered,]
  toInter_proportion_median <- means[intered,]
  toIner_proportion_change <- fold_changes[intered,]
  # time to make the formula
  values_with_preferences <- function(gene) {
    #print(gene)
    
    scaled_pref <- toInter_OR[gene,]
    gene_DE <- toInter_InGene[gene,]
    prop <- toInter_proportion_median[gene,]
    prop_fc <- toIner_proportion_change[gene, ]
    up <- gene_DE$log2fc > 0
    sign <- -1
    #print(up)
    if(up == T) {
      prop_fc <- 1/prop_fc
      sign <- 1
    } 
    val <- scaled_pref * prop * prop_fc * sign*2^abs(gene_DE$log2fc)
    return(val)
  }
  vals_out <- lapply(toInter_InGene$gene_name, values_with_preferences)
  vals_out_mat <- do.call("rbind", vals_out)
  rownames(vals_out_mat) <- toInter_InGene$gene_name
  
  # These are the genes that don't fit
  DE_out <- DEGs[genesOut,]
  all_gene_proportions <- proportion_pull(proportions)
  values_without_preferences <- function(gene) {
    scaled_pref <- 1
    gene_DE <- DEGs[gene,]
    prop <- all_gene_proportions$Mean
    prop_fc <- all_gene_proportions$FC
    up <- gene_DE$log2fc > 0
    if(up == T) {
      prop_fc <- 1/prop_fc
    } 
    val <- scaled_pref * prop * prop_fc * gene_DE$log2fc
    return(val)
  }
  vals_out_notpreffered <- lapply(genesOut, values_without_preferences)
  vals_out_notpreffered_mat <- do.call("rbind", vals_out_notpreffered)
  if(length(genesOut) == 0) {
    print("all genes are in signature matrix")
    return(vals_out_mat)
  }
  rownames(vals_out_notpreffered_mat) <- genesOut
  colnames(vals_out_notpreffered_mat) <- colnames(wilcox_or)
  colnames(vals_out_mat) <- colnames(wilcox_or)
  all_reordered <- rbind(vals_out_notpreffered_mat, vals_out_mat)
  comp <- plot_names
  myheatcol <- colorRampPalette(c("lightblue", "white", "orange"))(256)
  pdf(paste0(comp,"_heatmap.pdf"))
  cex = 0.4
  heatmap.2(data.matrix(vals_out_mat), Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
  dev.off()
  
  if(print_plots == T) {
    boxplot_values <- function(cmeaned_stacked, names) {
      all_stack <- c()
      for(i in 1:ncol(cmeaned_stacked)) {
        cm_stack <- cmeaned_stacked[,i]
        top3 <- head(sort(cm_stack), 3)
        bottom3 <- tail(sort(cm_stack), 3)
        empty <- rep("", length(cm_stack))
        names(empty) <- names(cm_stack)
        empty[names(top3)] <- names(top3)
        empty[names(bottom3)] <- names(bottom3)
        y <- cbind(names(cm_stack), unname(cm_stack), colnames(cmeaned_stacked)[i], unname(empty))
        all_stack <- rbind(all_stack, y)
      }
      all_stack <- as.data.frame(all_stack)
      colnames(all_stack) <- c("gene", "proportion", "cell_type", "label")
      all_stack$gene <- tochr(all_stack$gene)
      all_stack$proportion <- toNum(all_stack$proportion)
      all_stack$cell_type <- tochr(all_stack$cell_type)
      all_stack$label <- tochr(all_stack$label)
      library(ggplot2)
      pdf(paste0("deconvolute_generemove_quantseq_",names,".pdf"))
      g <-ggplot(all_stack, aes(factor(cell_type), proportion)) + geom_boxplot()  + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
      plot(g)
      dev.off()
      
      for(i in unique(all_stack$cell_type)) {
        cm_one <- all_stack[all_stack$cell_type == i,]
        pdf(paste0("deconvolute_generemov_quantseq_", i,"_",names,".pdf"))
        g <-ggplot(cm_one, aes(factor(cell_type), proportion)) + geom_boxplot() +geom_text(aes(label=label),hjust=-0.2) + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12))
        plot(g)
        dev.off()
        print(i)
      }
      return("Done!")
    }
    
    boxplot_values(cmeaned_stacked, paste0(plot_names,"_average"))
    boxplot_values(means, paste0(plot_names,"_median"))
    boxplot_values(fold_changes, paste0(plot_names,"_foldchange"))
  }
  
  return(all_reordered)
  
  
}

pathway_enrich <- function(i, name) {
  GL <- names(sort(abs(dm_test_nomt[,i][dm_test_nomt[,i] != 0]), decreasing = T))
  
  #unordered_noback <- gprofiler(GL, "hsapiens", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"))
  #unordered_back <- gprofiler(GL, "hsapiens", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
  ordered_back <- gprofiler(GL, "hsapiens", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3, hier_filtering = "moderate")
  l <- list( ordered_back = ordered_back, gene_list = GL)
  #names(l) <- name
  return(l)
}
pathway_enrich_TF <- function(i, name) {
  GL <- names(sort(abs(dm_test_nomt[,i][dm_test_nomt[,i] != 0]), decreasing = T))
  
  #unordered_noback <- gprofiler(GL, "hsapiens", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"))
  #ordered_noback <- gprofiler(GL, "hsapiens", ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),  min_isect_size = 3)
  #unordered_back <- gprofiler(GL, "hsapiens", ordered_query = F, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = rownames(norm_counts_i), min_isect_size = 3,  hier_filtering = "moderate")
  
  ordered_back <- gprofiler(GL, "hsapiens", ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = rownames(norm_counts_i), min_isect_size = 3, hier_filtering = "moderate")
  l <- list(  ordered_back = ordered_back, gene_list = GL)
  #names(l) <- name
  return(l)
}
