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
  
  norm_counts_i <- read.table(count_file, header = T, as.is = T, sep = "\t")
  load(odds_ratio)
  DEGs <- read.table(DEG_list, header = F, as.is = T, sep = "\t")
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
  all_genes_in <- DeconRNASeq(as.data.frame(norm_counts_i), wilcox_or)
  
  proportions <- all_genes_in$out.all
  rownames(proportions) <- colnames(norm_counts_i)
  
  propmeans <- colMeans(proportions)
  proportions <- proportions[,colSums(proportions) > 0]
  print("your bulk data contains the following cell types")
  print(colnames(proportions), quote = F)
  wilcox_or <- wilcox_or[,colnames(proportions)]
  wilcox_or_df <- as.data.frame(wilcox_or)
  bulk_in <- norm_counts_i
  genesIn <- intersect(DEGs$gene_name, rownames(wilcox_or))
  genesOut <- DEGs$gene_name[!(DEGs$gene_name %in% rownames(wilcox_or))]
  deconvolute_gene_removed <- function(x, bulk = bulk_in, signature = wilcox_or) {
    bulk_rem <- bulk[-which(rownames(bulk) == x),]
    signature_rem <- signature[-which(rownames(bulk) == x),]
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
    #print(up)
    if(up == T) {
      prop_fc <- 1/prop_fc
    } 
    val <- scaled_pref * prop * prop_fc * gene_DE$log2fc
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
  rownames(vals_out_notpreffered_mat) <- genesOut
  colnames(vals_out_notpreffered_mat) <- colnames(wilcox_or)
  colnames(vals_out_mat) <- colnames(wilcox_or)
  all_reordered <- rbind(vals_out_notpreffered_mat, vals_out_mat)
  comp <- 'reranked_pls'
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


library(DeconRNASeq)
library(robustbase)
library(gplots)
library(gProfileR)
count_file <- "C:/Users/Dustin Sokolowski/Documents/normalized_counts_symbol.txt"
odds_ratio <- "C:/Users/Dustin Sokolowski/Desktop/scMappR/heatmap/odd_ratio/merged_heatmap_or_SRA694785_Preoptic region of the hypothalamus.RData"
DEG_list <- "C:/Users/Dustin Sokolowski/Documents/quant_degs.tsv"
theSpecies <- -9
case_grep <- "37F"
control_grep <- "12F"
max_proportion_change <- -9
print_plots <- T
plot_names <- "test"
DEGs <- read.table(DEG_list, header = F, as.is = T, sep = "\t")
colnames(DEGs) <- c("gene_name", "padj", "log2fc")

norm_counts_i <- read.table(count_file, header = T, as.is = T,sep="\t")
tst <- deconvolute_and_contextualize(count_file, odds_ratio, DEG_list, case_grep = case_grep, control_grep = control_grep )
B <- Sys.time()
dm_test <- data.matrix(tst)

dm_test1 <- abs(dm_test)
dm_test2 <- dm_test1/min(dm_test1[dm_test1 != 0])
dm_test2[dm_test2 == 0] <- 1
dm_fin <- 1/dm_test2
dm_fin_nomt <- dm_fin[-grep("mt-",rownames(dm_fin)) ,]
#AP <- ActivePathways(dm_fin_nomt, "C:/Users/Dustin Sokolowski/Desktop/gmt/Mouse_GOBP_AllPathways_no_GO_iea_July_05_2019_symbol.gmt")

dm_test_nomt <- dm_test[-grep("mt-",rownames(dm_fin)),]
colnames(dm_test_nomt)
pathway_enrich <- function(i, name) {
  GL <- names(sort(abs(dm_test_nomt[,i][dm_test_nomt[,i] != 0]), decreasing = T))
  
  #unordered_noback <- gprofiler(GL, "mmusculus", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"))
  unordered_back <- gprofiler(GL, "mmusculus", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3)
  ordered_back <- gprofiler(GL, "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3)
  l <- list( unordered_back = unordered_back, ordered_back = ordered_back, gene_list = GL)
  #names(l) <- name
  return(l)
}
ordered_back <- gprofiler(geneList$Neuron, "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"), correction_method = "fdr", min_isect_size = 3)

pathway_enrich_TF <- function(i, name) {
  GL <- names(sort(abs(dm_test_nomt[,i][dm_test_nomt[,i] != 0]), decreasing = T))
  
  #unordered_noback <- gprofiler(GL, "mmusculus", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"))
  #ordered_noback <- gprofiler(GL, "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),  min_isect_size = 3)
  unordered_back <- gprofiler(GL, "mmusculus", ordered_query = F, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), min_isect_size = 3)
  ordered_back <- gprofiler(GL, "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = rownames(norm_counts_i), min_isect_size = 3)
  l <- list( unordered_back = unordered_back, ordered_back = ordered_back, gene_list = GL)
  #names(l) <- name
  return(l)
}
rownames(DEGs) <- DEGs$gene_name
DEGs_nomt <- DEGs[-grep("mt-",rownames(DEGs)),]
ordered_back_all <- gprofiler(DEGs_nomt$gene_name[order(DEGs_nomt$padj)], "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3)
ordered_back_nobg_all <- gprofiler(DEGs_nomt$gene_name[order(DEGs_nomt$padj)], "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"), correction_method = "fdr", min_isect_size = 3)
ordered_back_all_tf <- gprofiler(DEGs_nomt$gene_name[order(DEGs_nomt$padj)], "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3)

geneList <- list()
BP_endo <- list()
for(i in 1:ncol(tsts)){
endo <- pathway_enrich(i, "endo")
BP_endo[[i]] <- endo
geneList[[i]] <- endo$gene_list
print(colnames(dm_test_nomt)[i])
if(nrow(endo$ordered_back) > 0) {
  print(endo$ordered_back)
}
}


TF_endo <- list()
for(i in 1:ncol(tst)){
  endo <- pathway_enrich_TF(i, "endo")
  TF_endo[[i]] <- endo
  print(colnames(dm_test_nomt)[i])
  if(nrow(endo$ordered_back) > 0) {
    print(endo$ordered_back)
  }
}
C <- Sys.time()



# ### Ignore everything after this I'd say. I's just making boxplots of the p-vlaues andd is dumb.

plotBP <- function(ordered_back_all, top_bp = 10) {
ordered_back_all$term.name <- tolower(ordered_back_all$term.name)
ordered_back_all$log10 <- -1*log10(ordered_back_all$p.value)
if(nrow(ordered_back_all) > top_bp) {
  ordered_back_all <- ordered_back_all[1:top_bp,]
}
g <- ggplot(ordered_back_all, aes(x = reorder(term.name, log10), y = log10)) + geom_bar(stat = "identity", fill = "turquoise") + coord_flip() +  labs(y = "-log10(Padj)", x = "Gene Ontology") 
y <- g + theme(axis.text.x = element_text(face=NULL, color="black", 
                                    size=12, angle=35),
         axis.text.y = element_text(face=NULL, color="black", 
                                    size=12, angle=35), 
         axis.title=element_text(size=16, color = "black"))
print(y)
return(ordered_back_all)
}
make_TF_barplot <- function(ordered_back_all_tf, top_tf = 5) {
take1 <- function(x) return(x[1])
sp <- strsplit(ordered_back_all_tf$term.name, ";")
tfs <- unlist(lapply(sp, take1))
tfs <- gsub("Factor:","",gsub("-","", tochr(tfs)))
ordered_back_all_tf$tf <- tochr(tfs)
nodup <- ordered_back_all_tf[!duplicated(tochr(tfs)),]
ndup_1_10 <- nodup[order(nodup$p.value),]
if(nrow(ndup_1_10) > top_tf) {
ndup_1_10 <- ndup[1:top_tf,]
}
ndup_1_10$log10 <- -1*log10(ndup_1_10$p.value)
g <- ggplot(ndup_1_10, aes(x = reorder(tf, log10), y = log10)) + geom_bar(stat = "identity", fill = "mediumpurple") + coord_flip() +  labs(y = "-log10(Padj)", x = "TF Motif") 
y <- g + theme(axis.text.x = element_text(face=NULL, color="black", 
                                     size=12, angle=35),
          axis.text.y = element_text(face=NULL, color="black", 
                                     size=12, angle=35), 
          axis.title=element_text(size=16, color = "black"))
print(y)
return(ndup_1_10)
}
default_tf <- make_TF_barplot(ordered_back_all_tf)

CT <- c("Interneuron", "Neuron", "Oligo1", "endocrine/interneuron", "Oligo precursor",
        "Microglia", "Astro and Endocrine", "Neuroendocrine",
        "Oligo precursor1", "Schwann1","Endothelial","Mural","Oligo2",  "Astrocyte1")



names(BP_endo) <- names(TF_endo) <- names(geneList) <- CT


convert_em <- function(x, name) {
neuro_tf <- cbind(x$term.name, x$term.name, x$p.value/100, x$p.value/100, 1, x$intersection)
colnames(neuro_tf) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
write.table(neuro_tf, file = name, quote = F, row.names = F, col.names = T, sep = "\t")
return(neuro_tf)
}

convert_em <- function(x, name) {
  neuro_tf <- cbind(x$tf, x$tf, 1e-10, 1e-10, 1, x$intersection)
  colnames(neuro_tf) <- c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")
  write.table(neuro_tf, file = name, quote = F, row.names = F, col.names = T, sep = "\t")
  return(neuro_tf)
}



neuroendocrine_TF <- make_TF_barplot(TF_endo$Neuroendocrine$ordered_back, 10)
convert_em(neuroendocrine_TF, "neuroendocrine_tf.tsv")

oligo_TF <- make_TF_barplot(TF_endo$Oligo2$ordered_back, 10)
convert_em(oligo_TF, "oligo_tf.tsv")


neuron_TF <- make_TF_barplot(TF_endo$Neuron$ordered_back, 10)
convert_em(neuron_TF, "neuron_tf.tsv")


inter_endocrine_TF <- make_TF_barplot(TF_endo$`endocrine/interneuron`$ordered_back)
convert_em(inter_endocrine_TF, "inter_endocrine_tf.tsv")

neuroendocrine_BP <- plotBP(BP_endo$Neuroendocrine$ordered_back)
convert_em(neuroendocrine_BP, "neuroendocrine_BP.tsv")

Oligo_BP <- plotBP(BP_endo$Oligo2$ordered_back)
convert_em(Oligo_BP, "Oligo_BP.tsv")

Neuron_BP <- plotBP(ordered_back[1,])

convert_em(Neuron_BP, "Neuron_BP.tsv")

inter_endocrine_BP <- plotBP(BP_endo$`endocrine/interneuron`$ordered_back)
convert_em(inter_endocrine_BP, "inter_endocrine_BP.tsv")



