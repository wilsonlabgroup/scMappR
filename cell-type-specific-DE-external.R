library(DeconRNASeq)
library(robustbase)
library(gplots)
library(gProfileR)
load("data_for_validation.Rdata")
source("../cell-type-specific-DE-internal.R")
load("wilson_lab_examples/pituitary_background_preprocessed.Rdata")
counts <- readRDS("cadia_DE/pit_utr_2019__RUV_k2_set1_2019-07-03_HUAYUN.rds")
preprocessed_background[preprocessed_background < 0] = 0
# Load in the required files and get values for each cell type
count_file <- "C:/Users/Dustin Sokolowski/Documents/normalized_counts_symbol.txt"
odds_ratio <- preprocessed_background
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


# cell-type specific pathway analysis

B <- Sys.time()

# screwing with the output to make it ready for pathway analysis
dm_test <- data.matrix(tst)
dm_test_nomt <- dm_test[-grep("mt-",rownames(dm_test)),]

# pathway analysis -- no CT specificity
rownames(DEGs) <- DEGs$gene_name
DEGs_nomt <- DEGs[-grep("mt-",rownames(DEGs)),]
ordered_back_all <- gprofiler(DEGs_nomt$gene_name[order(DEGs_nomt$padj)], "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3)
ordered_back_nobg_all <- gprofiler(DEGs_nomt$gene_name[order(DEGs_nomt$padj)], "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 2000, src_filter = c("GO:BP", "REAC", "KEGG"), correction_method = "fdr", min_isect_size = 3)
ordered_back_all_tf <- gprofiler(DEGs_nomt$gene_name[order(DEGs_nomt$padj)], "mmusculus", ordered_query = T, min_set_size = 3, max_set_size = 5000, src_filter = c("TF"),custom_bg = rownames(norm_counts_i), correction_method = "fdr", min_isect_size = 3)
##

# GO:BP analysis of each cell-type
geneList <- list()
BP_endo <- list()
for(i in 1:ncol(tst)){
endo <- pathway_enrich(i, "endo")
BP_endo[[i]] <- endo
geneList[[i]] <- endo$gene_list
print(colnames(dm_test_nomt)[i])
if(nrow(endo$ordered_back) > 0) {
  print(endo$ordered_back)
}
}

# TF enrichment analysis of each cell type
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



