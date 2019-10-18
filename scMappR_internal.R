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

character_correlate <- function(x) {
  if(class(x) == "factor") return(as.numeric(x))
  if(class(x) == "character") return(as.numeric(factor(x)))
}



single_gene_preferences <- function(hg_short, hg_full, study_name, outDir = output_directory) {

fP <- hg_full$preferences
sP <- hg_short$preferences
in_length <- length(hg_short$genesIn)
inmarker_length <- length(hg_full$genesIn)
fP <- hg_full$preferences
sP <- hg_short$preferences
in_length <- length(hg_short$genesIn)

pref <- c()
for(cell in 1:length(sP)) {
theGenes <- unique(sP[[cell]])
theGenes_p <- paste0(theGenes,collapse = ",")
aa <- length(unique(sP[[cell]]))
ab <- in_length - aa
ba <- length(unique(fP[[cell]]))
bb <-  inmarker_length - ba
m <- matrix(c(aa, ba, ab, bb), nrow = 2)
ftest <- fisher.test(m)
p <- ftest$p.value
OR <- ftest$estimate
nSP <- names(sP)[cell]
P <- c(nSP, p, OR, theGenes_p)
pref <- rbind(pref, P)
}
colnames(pref) <- c("cell_type", "p_val", "Odds_Ratio", "genes")
pref <- as.data.frame(pref)
write.table(pref, file = paste0(outDir, "/",study_name, "cell_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
pref <- read.table(file = paste0(outDir, "/",study_name, "cell_preferences.tsv"), as.is=T,header=T, sep = "\t")
pref$pFDR <- p.adjust(pref$p_val, "fdr")
write.table(pref, file = paste0(outDir, "/",study_name, "cell_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")

return(pref)
}
extract_genes_cell <- function(geneHeat, cellTypes = "ALL", pVal = 0.01, isMax = F) {
#geneHeat: the heatmap of ranks from your scRNA-seq dataset with your genes subsetted
#celltyes: The cell-types that you're interested in extracting. They need to be colnames (not-case sensitive)
# pVal, how associated a gene is with a particualr cell type to include in your list
# isMax true or false, if you want to sort genes into the CT marker that represents theme best


colnames(geneHeat) <- toupper(colnames(geneHeat))
cellTypes <- toupper(cellTypes)

geneHeat <- geneHeat[rowSums(geneHeat) > 1,] # there is some preference at all
if(class(geneHeat) == "numeric" | class(geneHeat) == "character" ) {
  geneHeat <- data.frame(t(geneHeat))
}  

genes_extracted <- list()
if(cellTypes == "ALL") {
  cellTypes <- colnames(geneHeat)
}

if(isMax == TRUE) { 
  # this way, you take the cell type that marks each gene the best.
  # no gene is double counted into multiple cell types.
  cn <- colnames(geneHeat)[max.col(geneHeat,"first")] # the cell type with the most signficant
  names(cn) <- rownames(geneHeat)
  for(i in 1:length(colnames(geneHeat))) {
    genes_extracted[[i]] <- names(cn)[which(cn == colnames(geneHeat)[i])]
    names(genes_extracted)[i] <- colnames(geneHeat)[i]
  }
  return(genes_extracted) 
} else {
  # a gene can be counted into multiple cell types and you're just saying that a 
  # gene is in it if it's p-value passes a certain threshold -- This is probably
  # what fits the premise of scMappR slightly better
  for(i in 1:length(colnames(geneHeat))) {
    genes_extracted[[i]] <- rownames(geneHeat)[geneHeat[,i] > -1*log10(pVal)]
    names(genes_extracted)[i] <- names(genes_extracted)[i] <- colnames(geneHeat)[i]
    
  }
  return(genes_extracted)
}
}

heatmap_generation <- function(genesIn, comp, cex = genecex, cellTypes = "ALL", pVal = 0.01, isMax =F,  isBackground = F,reference = "C:/Users/Dustin Sokolowski/Desktop/romanov_wilcoxon_test_2.RData", custom = FALSE, which_species = -9) {
#Hi Mari! There are some parameters for this guy :)
# genesIn = a list of gene symbols (all caps) to have their cell type enrichment
# comp = the name of the comparison
# cex = the size of the genes in the column label for the heatmap... I know
# The parameters below are passed into extract_genes_cell
# cellTypes = which cell types you care about extracting (colnames)
# For romanov, it's as follows: "OLIGOS", "ASTROCYTES", "EPENDYMAL", "MICROGLIA", "VSM", "ENDOTHELIAL", "ADCYAP1", "AVP", "DOPAMINE", "GABA", "OTHER", "VGLUT", "CIRCADIAN"  
# pVale = the level of association a gene has within a cell type
# isMax = true/false, just sort the gene into the cell type it's most strongly associated with
library(gplots)
if(custom == TRUE) {
  wilcoxon_map_rnk_t <- reference  
  
} else {
  load(reference) # load the summary stats
  wilcoxon_map_rnk_t <- wilcoxon_map_rnk_t[!duplicated(rownames(wilcoxon_map_rnk_t)),]
  if(length(grep("-", rownames(wilcoxon_map_rnk_t))) > 50) { 
    
    RN_2 <- get_gene_symbol(wilcoxon_map_rnk_t)
    rownames(wilcoxon_map_rnk_t) <- RN_2$rowname
  }
  if(which_species != RN_2$species ) {
    ### TERMPORARY FIX
    rownames(wilcoxon_map_rnk_t) <- toupper(rownames(wilcoxon_map_rnk_t))
    genesIn <- toupper(genesIn)
  }
  
}
if(any(duplicated(colnames(wilcoxon_map_rnk_t)))) {
  print("cell-types not uniquely named, appending tag to make all cell-types unique")
  colnames(wilcoxon_map_rnk_t) <- paste0(colnames(wilcoxon_map_rnk_t),"_tag",1:ncol(wilcoxon_map_rnk_t))
  
}
wilcoxon_map_rnk_t <- wilcoxon_map_rnk_t[!duplicated(rownames(wilcoxon_map_rnk_t)),]

genesInter <- intersect(genesIn, rownames(wilcoxon_map_rnk_t)) # get genes with some DE
whichGenesInter <- which(rownames(wilcoxon_map_rnk_t) %in% genesInter)
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
  pdf(paste0(comp,"_barplot.pdf"))
  barplot(wilcoxon_map_rnk_t[whichGenesInter,], las = 2, main = rownames(wilcoxon_map_rnk_t)[whichGenesInter])
  dev.off()
  return("No_downstream_analysis")
}
if(length(whichGenesInter) > 1) { # if > 1 genes are
  print("at least one input gene is preferentially expressed")
  # make the heatmap
  myheatcol <- colorRampPalette(c("lightblue", "white", "orange"))(256)
  pdf(paste0(comp,"_heatmap.pdf"))
  heatmap.2(wilcoxon_map_rnk_t[whichGenesInter,], Rowv = T, dendrogram = "column", col = myheatcol, scale = "row", trace = "none", margins = c(7,7),cexRow = cex, cexCol = 0.3 )
  dev.off()
  geneHeat <- wilcoxon_map_rnk_t[whichGenesInter,]
  # write the table that made the heatmap (rank of every gene in every cell)
  #write.table(geneHeat, file = paste0(comp,"_gene_CT_rank_mat.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
  # get lists of cell-type preferences
  #print(head(comp))
  preferences <- extract_genes_cell(geneHeat, cellTypes = cellTypes, pVal = pVal, isMax = isMax)
  save(preferences, file = paste0(comp,"_preferences.RData"))
} 


return(list(genesIn = genesInter, genesNoIn = genes_noInter, geneHeat = geneHeat, preferences = preferences))
# genesIn = genes with some preference, genesNoIn = genes with no preference 
# geneHeat = matrix of preferences and p-values
# preferences = genes mapped to their CT preference
}


coEnrich <- function(sig, gene_list_heatmap, background_heatmap, study_name, outDir = output_directory) {
if(nrow(sig) > 5) {
  sig <- sig[1:5,]
}
l <- nrow(sig)
multi_comps <- c()
for(y in 2:l) {
  if(y < l) {
    comps <- combn(sig$cell_type, y)
    co_up <- function(x) return(length(x[x>=1])==y)
    print(y)
    for(j in 1:ncol(comps)) {
      thecomps <- comps[,j]
      geneList_comb <- gene_list_heatmap$geneHeat
      colnames(geneList_comb) <- toupper(colnames(geneList_comb))
      geneList_comb1 <- geneList_comb[,thecomps]
      background_comb <- background_heatmap$geneHeat
      colnames(background_comb) <- toupper(colnames(background_comb))
      background_comb1 <- background_comb[-which(rownames(background_comb) %in% rownames(geneList_comb)), thecomps]
      bin_aa <-apply(geneList_comb1,1,co_up)
      name_in <- names(bin_aa)[bin_aa == T]
      aa <- sum(bin_aa)
      ab <- nrow(geneList_comb1) - aa
      ba <- sum(apply(background_comb1,1,co_up))
      bb <- nrow(background_comb1) - ba
      m <- matrix(c(aa, ba, ab, bb), nrow = 2)
      fisherTest <- fisher.test(m)
      OR <- fisherTest$estimate
      p <- fisherTest$p.value
      name <- paste0(thecomps,collapse=":" )
      row <- c(name,p,OR, paste0(name_in, collapse = ","))
      multi_comps <- rbind(multi_comps, row)
    }
  }
  if(y == l ) {
    co_up <- function(x) return(length(x[x>=1])==y)
    
    thecomps <- sig$cell_type
    geneList_comb <- gene_list_heatmap$geneHeat
    colnames(geneList_comb) <- toupper(colnames(geneList_comb))
    geneList_comb1 <- geneList_comb[,thecomps]
    background_comb <- background_heatmap$geneHeat
    colnames(background_comb) <- toupper(colnames(background_comb))
    background_comb1 <- background_comb[-which(rownames(background_comb) %in% rownames(geneList_comb)), thecomps]
    bin_aa <-apply(geneList_comb1,1,co_up)
    name_in <- names(bin_aa)[bin_aa == T]
    aa <- sum(bin_aa)
    ab <- nrow(geneList_comb1) - aa
    ba <- sum(apply(background_comb1,1,co_up))
    bb <- nrow(background_comb1) - ba
    m <- matrix(c(aa, ba, ab, bb), nrow = 2)
    fisherTest <- fisher.test(m)
    OR <- fisherTest$estimate
    p <- fisherTest$p.value
    name <- paste0(thecomps,collapse=":" )
    row <- c(name,p,OR, paste0(name_in, collapse = ","))
    multi_comps <- rbind(multi_comps, row)
  }
}
colnames(multi_comps) <- c("cell_types", "p_val", "OR", "genes")
write.table(multi_comps, file = paste0(outDir, "/",study_name, "cell_co_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
multi_comps <- read.table(file = paste0(outDir, "/",study_name, "cell_co_preferences.tsv"), as.is=T,header=T, sep = "\t")
multi_comps$pFDR <- p.adjust(multi_comps$p_val, "fdr")
write.table(multi_comps, file = paste0(outDir, "/",study_name, "cell_co_preferences.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
return(multi_comps)
}


get_tissues <- function() {
  hm <- list.files(pattern = "*.RData", path = "heatmap/p_val/")
  hm1 <- substr(hm, 1, nchar(hm)-6)
  hm2 <- substr(hm1, 28, nchar(hm1))
  uni <- sort(unique(hm2))
  return(uni)
}

tissue_scMappR_internal <- function(gene_list,species, output_directory, tissue, cluster = "Pval", genecex = 0.01) {
  outDir <- output_directory
  #scMappR_default_dataset <- function(genesIn, comp, cex = 0.05, cellTypes = "ALL", pVal = 0.01, isMax =F,  isBackground = F, species = "mouse")
  hm <- list.files(pattern = "*.RData", path = "heatmap/p_val/")
  hm_tissue <- grep(toupper(tissue), toupper(hm))
  input_studies <- hm[hm_tissue]
  study_names <- substr(input_studies, 18, nchar(input_studies) - 6)
  dir.create(outDir)
  single_cell_studies <- list()
  for(i in 1:length(study_names)) {
    study_ref <- paste0("heatmap/p_val/",input_studies[i])
    load(paste0("heatmap/p_val/",input_studies[i]))
    sym <- get_gene_symbol(wilcoxon_map_rnk_t)
    
    
    background_genes <- sym$rowname
    theSpecies <- sym$species
    #print(background_genes)
    #test_genes <- sample(background_genes, 100, F)
    background_heatmap <- heatmap_generation(background_genes, comp = paste0(outDir, "/", study_names[i],"_background"), reference = study_ref, isBackground = TRUE, cex = genecex, which_species = species)  
    gene_list_heatmap <- heatmap_generation(gene_list, comp = paste0(outDir, "/", study_names[i],"_genelist"), reference = study_ref, cex = genecex, which_species = species)
    #print(class(gene_list_heatmap))
    if(class(gene_list_heatmap) == "character") {
      print("not enough genes were present to do downsteam analysis in: ")
      print(study_names[i])
      single_cell_studies[[i]] <- gene_list_heatmap
      next
    }
    singleCTpreferences <- single_gene_preferences(gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory)
    sig <- singleCTpreferences[singleCTpreferences$pFDR < 0.05,]
    sig <- sig[sig$Odds_Ratio > 1,]
    
    if(nrow(sig) <2 ) {
      print("co-enrichment cannot be measured as one or fewer CTs are enriched")
      coCTpreferences <- "co-enrichment cannot be measured as one or fewer CTs are enriched"
      output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
      single_cell_studies[[i]] <- output 
      next
    }
    sig <- sig[order(sig$pFDR),]
    coCTpreferences <- coEnrich(sig, gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory)
    output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
    
    single_cell_studies[[i]] <- output 
  }
  names(single_cell_studies) <- study_names
  return(single_cell_studies)
}



#
################# Now when you have a custom set
#
tissue_scMappR_custom <- function(gene_list, background ,output_directory, study_name) {
  outDir <- paste0(output_directory)
  #scMappR_default_dataset <- function(genesIn, comp, cex = 0.05, cellTypes = "ALL", pVal = 0.01, isMax =F,  isBackground = F, species = "mouse")

  study_names <- study_name
  dir.create(outDir)
  single_cell_studies <- list()
  for(i in 1:length(study_names)) {
    #test_genes <- sample(background_genes, 100, F)
    background_genes <- rownames(background)
    background_heatmap <- heatmap_generation(background_genes, comp = paste0(outDir, "/", study_names[i],"_background"), reference = background, isBackground = TRUE, custom = TRUE)  
    gene_list_heatmap <- heatmap_generation(gene_list, comp = paste0(outDir, "/", study_names[i],"_genelist"), reference = background, custom = TRUE )
    singleCTpreferences <- single_gene_preferences(gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory)
    sig <- singleCTpreferences[singleCTpreferences$pFDR < 0.05,]
    sig <- sig[sig$Odds_Ratio > 1,]
    
    if(nrow(sig) <2 ) {
      print("co-enrichment cannot be measured as one or fewer CTs are enriched")
      coCTpreferences <- "co-enrichment cannot be measured as one or fewer CTs are enriched"
      output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
      single_cell_studies[[i]] <- output 
      next
    }
    sig <- sig[order(sig$pFDR),]
    coCTpreferences <- coEnrich(sig, gene_list_heatmap, background_heatmap, study_names[i], outDir = output_directory)
    output <- list(background_heatmap = background_heatmap, gene_list_heatmap = gene_list_heatmap, single_celltype_preferences = singleCTpreferences, group_celtype_preference = coCTpreferences)
    
    single_cell_studies[[i]] <- output 
  }
  names(single_cell_studies) <- study_names
  return(single_cell_studies)
}

topgenes_extract <- function(generes,  padj = 0.05, FC = 1.5) {
  # This function takes a generes object and extracts the top (up to) 30 CT markers
  # it then returns it as a list
  topGenes <- list()
  for(i in 1:length(generes)) {
    genes <- rownames(generes[[i]])[generes[[i]]$avg_logFC > log2(FC) & generes[[i]]$p_val_adj < padj]
    if(length(genes) > 30) {
      genes <- genes[1:30]
    }
    topGenes[[i]] <- genes
  }
  names(topGenes) <- names(generes)
  return(topGenes)  
}

cellmarker_enrich <- function(toTest_n, p_thresh, gmt = "cellmarker_list.Rdata", fixed_length = 13000) {
  tochr <- function(x) {
    if(class(x) == "charater") return(x)
    if(class(x) == "factor") return(as.character(levels(x))[x])
    if(class(x) == "numeric") return(as.character(x))
  } 
  
  toNum <- function(x) {
    if(class(x) == "character") return(as.numeric(x))
    if(class(x) == "factor") return(as.numeric(levels(x))[x])
    if(class(x) == "numeric") return(x)
  } 
  load(gmt)
  theRows <- c()
  for(i in 1:length(gmt)) {
    
    test <-gmt[[i]]
    if(length(test) < 1 | length(test) > 3000) {
      next
    }
    ints <- intersect(toTest_n, test)
    tbt = matrix(c(length(ints), length(toTest_n)-length(ints),  length(test) ,fixed_length-length(test)), nrow=2)
    #phyper(length(ints), length(test),  length(toTest_n), 58037-length(test))
    rownames(tbt) <- c("Exposure_yes", "Exposure_no")
    colnames(tbt) <- c("outcome_yes", "outcome_no")
    f_out <- fisher.test(tbt, B = 1e9)
    p <- f_out$p.value
    OR <- f_out$estimate[[1]]
    if(p < 0.05 & OR < 1) {
      p = 1
    }
    term_size <- length(test)
    intersect_size <- length(ints)
    input_length <- length(toTest_n)
    proportion_success <- signif(intersect_size / term_size,2)
    interGenes <- paste0(ints, collapse = ",")
    theName <- names(gmt)[i]
    theRow <- c(theName, p, term_size, intersect_size, input_length, interGenes) 
    theRows <- rbind(theRows, theRow)
  }
  colnames(theRows) <- c("name","p", "term_size", "intersect_size", "input_length", "genes")
  df_theRows <- as.data.frame(theRows)
  if(class(df_theRows$p) == "factor") {
    df_theRows$p <- toNum(df_theRows$p)
  }
  if(class(df_theRows$p) == "character") {
   df_theRows$p <- as.numeric(df_theRows$p) 
  }
  fdr <- p.adjust(df_theRows$p, "fdr")
  bonf <- p.adjust(df_theRows$p, "bonferroni")
  df_theRows$fdr <- fdr
  df_theRows$bonf <- bonf
  if(class(df_theRows$intersect_size) == "factor") {
    df_theRows$intersect_size <- toNum(df_theRows$intersect_size)
  }
  if(class(df_theRows$intersect_size) == "character") {
    df_theRows$intersect_size <- as.numeric(df_theRows$intersect_size)
  }
  
  df_theRows_order <- df_theRows[order(df_theRows$p),]
  
  #rownames(df_theRows_order) <- tochr(df_theRows_order$name)
  df_theRows_order <- df_theRows_order[df_theRows_order$intersect_size > 2,]
  
  sig_studies <- df_theRows_order[df_theRows_order$p < p_thresh,]
  return(sig_studies)
}


human_mouse_ct_marker_enrich <- function(gene_lists, theSpecies = "human",naming_preference = -9) {
  topGenes <- gene_lists
  marker_sets <- list()
  
  tochr <- function(x) {
    # make sure that the vector is of character
    if(class(x) == "character") return(x)
    if(class(x) == "factor") return(as.character(levels(x))[x])
    if(class(x) == "numeric") return(as.character(x))
  } 
  
  toNum <- function(x) {
    # same deal but numeric
    if(class(x) == "character") return(as.numeric(x))
    if(class(x) == "factor") return(as.numeric(levels(x))[x])
    if(class(x) == "numeric") return(x)
  }   
  if(class(topGenes)== "character") {
    topGenes1 <- list()
    topGenes1[[1]] <- topGenes
    topGenes <- topGenes1
  }
  
  cellTypes <- c()
  
  for(i in 1:length(topGenes)) {
    #print(i)
    geneInput <- topGenes[[i]]
    if(theSpecies == "mouse") { # if it's a mouse
      
      outCellMarker <- cellmarker_enrich(topGenes[[i]], 0.05, "cell_type_info/cellmarker_list_mouse.Rdata", fixed_length = 15000)
      # get the cell type from 'cellmarkers' dataset
      if(class(outCellMarker$name) == "factor") {
        outCellMarker$name <- tochr(outCellMarker$name)
      }
      outPanglao <- cellmarker_enrich(topGenes[[i]], 0.05, "cell_type_info/mouse_marker_panglao.RData", fixed_length = 15000)
      # get the cell type from 'panglao' dataset
      if(class(outPanglao$name) == "factor") {
        outPanglao$name <- tochr(outPanglao$name)
      }
      outGOMouse <- cellmarker_enrich(topGenes[[i]], 0.05, "cell_type_info/mouse_GOBP.Rdata", fixed_length = 22000)
      # get the cell type from gene ontology
      if(class(outGOMouse$name) == "factor") {
        outGOMouse$name <- tochr(outGOMouse$name)
      }
      
      marker_list <- list(CellMarker = outCellMarker, Panglao = outPanglao, Ontology = outGOMouse)
      # concatenate all of the CT markers from the appropriate dataset
      if(naming_preference != -9) {
        load("cell_type_info/cell_preferences_categorized.Rdata")
        # if they have the prefered tissue or cell-type to pick from
        # prioritize those
        mypref <- cell_preference_final[[naming_preference]]
        outcell_preference <- which(outCellMarker$name %in% mypref)
        if(length(outcell_preference) > 0) {
          outCellMarker <- outCellMarker[outcell_preference,]
        }
        panglao_preference <- which(outPanglao$name %in% mypref)
        if(length(panglao_preference) > 0) {
          outPanglao <- outPanglao[panglao_preference,]
        }
        
      }
      marker_sets[[i]] <- marker_list 
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) == 0) {
        cellTypes[i] <- "unknown"
        # if none of that works out it's unknown
      }
      
      if(nrow(outGOMouse) == 0 & nrow(outGOMouse) == 0 & nrow(outGOMouse) > 0) {
        thenames <- sub( "\\%.*", "", outGOMouse$name)
        cellTypes[i] <- paste0("Ontology ", thenames[1])
        # define by top gene ontology if no other marker is present
        
      }
      
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) == 0) {
        thenames <- sub('.*\\:', '', outCellMarker$name)
        cellTypes[i] <- thenames[1]
        # define by cellmarker alone if none of the others are present
      }
      
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) > 0) {
        thenames <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- thenames[1]
        # define by panglao if it's the only dataset with a marker
      }
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) > 0) {
        thenamesCM <- sub('.*\\:', '', outCellMarker$name)
        thenamesP <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- paste0(thenamesCM[1], "_and_", thenamesP[1])
        # define by both if there was a signfiicant markerin both
      }
      
    }
    if(theSpecies == "human") {
      # this section is exactly the same as in mouse but generated for human
      # cell-marker lists
      outCellMarker <- cellmarker_enrich(topGenes[[i]], 0.05, "cell_type_info/cellmarker_list_human.Rdata", fixed_length = 15000)
      if(class(outCellMarker$name) == "factor") {
        outCellMarker$name <- tochr(outCellMarker$name)
      }
      outPanglao <- cellmarker_enrich(topGenes[[i]], 0.05, "cell_type_info/human_marker_panglao.RData", fixed_length = 15000)
      if(class(outPanglao$name) == "factor") {
        outPanglao$name <- tochr(outPanglao$name)
      }
      outGoHuman <- cellmarker_enrich(topGenes[[i]], 0.05, "cell_type_info/human_GOBP.Rdata", fixed_length = 22000)
      if(class(outGoHuman$name) == "factor") {
        outGoHuman$name <- tochr(outGoHuman$name)
      }
      if(naming_preference != -9) {
        load("cell_type_info/cell_preferences_categorized.Rdata")
        mypref <- cell_preference_final[[naming_preference]]
        outcell_preference <- which(outCellMarker$name %in% mypref)
        if(length(outcell_preference) > 0) {
          outCellMarker <- outCellMarker[outcell_preference,]
        }
        panglao_preference <- which(outPanglao$name %in% mypref)
        if(length(panglao_preference) > 0) {
          outPanglao <- outPanglao[panglao_preference,]
        }
        
      }
      marker_list <- list(CellMarker = outCellMarker, Panglao = outPanglao, Ontology = outGoHuman)
      marker_sets[[i]] <- marker_list 
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) == 0 & nrow(outGoHuman) == 0) {
        cellTypes[i] <- "unknown"
      }
      
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) == 0 & nrow(outGoHuman) > 0) {
        thenames <- sub( "\\%.*", "", outGoHuman$name)
        cellTypes[i] <- paste0("Ontology ", thenames[1])
        
      }
      
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) == 0) {
        thenames <- sub('.*\\:', '', outCellMarker$name)
        cellTypes[i] <- thenames[1]
      }
      
      if(nrow(outCellMarker) == 0 & nrow(outPanglao) > 0) {
        thenames <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- thenames[1]
      }
      if(nrow(outCellMarker) > 0 & nrow(outPanglao) > 0) {
        thenamesCM <- sub('.*\\:', '', outCellMarker$name)
        thenamesP <- gsub("_panglao","",outPanglao$name)
        cellTypes[i] <- paste0(thenamesCM[1], "_and_", thenamesP[1])
      }
      
    }
  }
  theNeurons <- grep("Neuron", cellTypes)
  # for neuronal cell-types take a secondary round of classification to see if it's a neuronal subtype
  if(length(theNeurons) > 0) {
    TGenes <- topGenes[theNeurons]
    #neuron_name <- names(TGenes)
    
    TG <- c()
    for(z in 1:length(TGenes)) {
      if(theSpecies == "mouse") {
        outNeuronalMouse <- cellmarker_enrich(TGenes[[z]], 1, "cell_type_info/mouse_neuronal_subtype.Rdata", fixed_length = 15000)
        
        if(nrow(outNeuronalMouse) == 0) {
          TG[z] <- "neuron_unclassified"
        } else {
          if(class(outNeuronalMouse$name) == "factor") {
            outNeuronalMouse$name <- tochr(outNeuronalMouse$name)
          }
          TG[z] <- unname(unlist(outNeuronalMouse$name[1]))
        }
        #print(outNeuronalMouse)
      }
      if(theSpecies == "human") {
        outNeuronalHuman <- cellmarker_enrich(TGenes[[z]], 1, "cell_type_info/human_neuronal_subtype.Rdata", fixed_length = 15000)
        if(nrow(outNeuronalHuman) == 0) {
          TG[z] <- "neuron_unclassified"
        } else {
          if(class(outNeuronalHuman$name) == "factor") {
            outNeuronalHuman$name <- tochr(outNeuronalHuman$name)
          }
          TG[z] <- unname(unlist(outNeuronalHuman$name[1]))
        }
        #print(outNeuronalMouse)
      }
    }
    cellTypes[theNeurons] <- TG
  }
  names(marker_sets) <- cellTypes
  return(list(cellTypes = cellTypes, marker_sets = marker_sets))
}



generes_to_heatmap <- function(generes,  species = "human", name = "test", naming_preference = -9) {
  topGenes <- topgenes_extract(generes) # take the top 30 genes
  cell_name <- human_mouse_ct_marker_enrich(topGenes, theSpecies = species) 
  # get the names of each of the cell types
  cell_name <- cell_name$cellTypes # attach the appropriate cell names
  genes <- c()
  for(i in 1:length(generes)) {
    genes <- c(genes,rownames(generes[[i]]))
  }
  # generate matrix with ranks or odds ratio
  # only genes that are preferentially expressed in at min one cell type
  # is removed
  genes_uni <- unique(genes)
  scmappr <- matrix(0, nrow = length(genes_uni), ncol = length(names(generes)))
  rownames(scmappr) <- genes_uni
  colnames(scmappr) <- names(generes)
  for(i in 1:length(generes)) {
    rnk <- -1*log10(generes[[i]]$p_val_adj) * sign(generes[[i]]$avg_logFC)
    names(rnk) <- rownames(generes[[i]])
    scmappr[names(rnk),i] <- unname(rnk)
  }
  wilcoxon_map_rnk_t <- scmappr
  wilcoxon_map_rnk_t[is.infinite(wilcoxon_map_rnk_t) & wilcoxon_map_rnk_t < 0] <- min(wilcoxon_map_rnk_t[is.finite(wilcoxon_map_rnk_t)])
  wilcoxon_map_rnk_t[is.infinite(wilcoxon_map_rnk_t) & wilcoxon_map_rnk_t > 0] <- max(wilcoxon_map_rnk_t[is.finite(wilcoxon_map_rnk_t)])
  
  scmappr <- matrix(0, nrow = length(genes_uni), ncol = length(names(generes)))
  rownames(scmappr) <- genes_uni
  colnames(scmappr) <- names(generes)
  for(i in 1:length(generes)) {
    rnk <- generes[[i]]$avg_logFC
    names(rnk) <- rownames(generes[[i]])
    scmappr[names(rnk),i] <- unname(rnk)
  }
  wilcoxon_map_rnk_or <- scmappr
  wilcoxon_map_rnk_or[is.infinite(wilcoxon_map_rnk_or) & wilcoxon_map_rnk_or < 0] <- min(wilcoxon_map_rnk_or[is.finite(wilcoxon_map_rnk_or)])
  wilcoxon_map_rnk_or[is.infinite(wilcoxon_map_rnk_or) & wilcoxon_map_rnk_or > 0] <- max(wilcoxon_map_rnk_or[is.finite(wilcoxon_map_rnk_or)])
  l <- list(pVal = wilcoxon_map_rnk_t, OR = wilcoxon_map_rnk_or, cellname = cell_name)
  return(l)
}

seurat_to_generes <- function(pbmc){
  # cell identity within seurat object must be named ident (seurat 2) or acive.ident (seurat 3)
  id <- try(pbmc@ident, silent = T)

  generes <- list()
  count <-1
  if(class(id) != "try-error") {
    for(i in sort(unique(pbmc@ident))) {
      #toUse <- as.character(pheno_df_o[i,])
      #toUse <- toUse[pbmc@cell.names]
      #toUse <- as.factor(toUse)
      #pbmc@ident <- toUse
      de_genes <- try(FindMarkers(pbmc, ident.1 = i, test.use = "wilcox"))
      if(class(de_genes) == "try-error" | class(de_genes) == "NULL" ) {
        
        next
      }
      
      generes[[count]] <- de_genes
      names(generes)[count] <- i
      count <- count+1
      print("DE")
      print(i)
    }
    
  }
  for(i in sort(unique(pbmc@active.ident))) {
    #toUse <- as.character(pheno_df_o[i,])
    #toUse <- toUse[pbmc@cell.names]
    #toUse <- as.factor(toUse)
    #pbmc@ident <- toUse
    de_genes <- try(FindMarkers(pbmc, ident.1 = i, test.use = "wilcox"))
    if(class(de_genes) == "try-error" | class(de_genes) == "NULL" ) {
      
      next
    }
    
    generes[[count]] <- de_genes
    names(generes)[count] <- i
    count <- count+1
    print("DE")
    print(i)
  }
  return(generes)
}

process_from_count <- function(countmat_list, name, theSpecies = -9, haveUmap = FALSE, saveALL = FALSE, panglao_set = F) {
  
  library("Seurat") 
  library("sctransform")
  SRA_in <- countmat_list
  
  shrt <- names(SRA_in)
  name <- name
  each_sra <- list()
  count <- 1
  for(f in 1:length(SRA_in)) {
    
    sm <- SRA_in[[f]]
    sm <- sm[!duplicated(rownames(sm)),]
    ####################################
    ####################################
    
    if(theSpecies == -9 | panglao_set == TRUE) {
      
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
      print(theSpecies)
      # make sure that mitochondrial genes are flagged
      RN = rownames(sm)
      RN_1 <- sub('(.*)[.](.*)','\\1',RN)
      if(theSpecies == "mouse") {
      RN_2 <- substr(RN_1, 1, nchar(RN_1) - 19)
      }
      if(theSpecies == "human") {
        RN_2 <- substr(RN_1, 1, nchar(RN_1) - 16)
      }
      rownames(sm) <- RN_2

    }
    RN_2 <- rownames(sm)
    if(theSpecies == "mouse") {
      num_MT <- grep("mt-", RN_2)
    }
    if(theSpecies == "human") {
      num_MT <- grep("MT-",RN_2)
    }
    
    # the code below will check to see if there are mt genes in the correct format.
    # If they are in the right format then continue on, 
    # otherwise adjust the gene names so that mito genes are detected
    mt.genes_m <- c("Tf", "Rnr1","Tv","Rnr2","Tl1","Nd1","Ti","Tq","Tm","Nd2","Tw","Ta","Tn","Tc","Ty","Co1","Ts1","Td","Co2","Tk","Atp8","Atp6","Co3","Tg","Nd3","Tr","Nd4l","Nd4","Ts2","Tl2","Nd5","Nd6","Te","Cytb","Tt","Tp")
    print(length(num_MT))
    if(length(num_MT) == 0 & theSpecies =="human") {
      # convert gene names
      mt.genes <- toupper(mt.genes_m)
      mito.genes <- which(RN_2 %in% mt.genes)
      RN_2[mito.genes] <- paste0("MT-", RN_2[mito.genes])
      
      mito.genes <- RN_2[mito.genes]
    }
    
    if(length(num_MT) == 0 & theSpecies =="mouse") {
      # convert gene names mouse
      mt.genes <- mt.genes_m
      mito.genes <- which(RN_2 %in% mt.genes)
      mito.genes <- RN_2[mito.genes]
      RN_2[mito.genes] <- paste0("mt-", RN_2[mito.genes])
      
    }
    rownames(sm) <- RN_2
    
    pbmc <- Seurat::CreateSeuratObject(sm, min.cells = 3, min.features = 0, project = shrt[count])
    if(theSpecies == "human") {
      pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt.adj")
    }
    if(theSpecies == "mouse") {
      pbmc <- PercentageFeatureSet(pbmc, pattern = "^mt-", col.name = "percent.mt.adj")
    }
    
    # TEMPORARY!
    #########
    #pbmc <- pbmc[,which(pbmc$percent.mt.adj < 50)]
    mean <- mean(pbmc$percent.mt.adj)
    x<- sd(pbmc$percent.mt.adj)
    toremove <- mean + (2*x)
    pbmc <- pbmc[,which(pbmc$percent.mt.adj < toremove)]
    #########
    #
    # Cannot process with scTransform because it doesn't allow for integration yet -- may be able to update
    pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt.adj", verbose = FALSE)
    # process with sctransformed (most updated)
    
    
    
    each_sra[[count]] <- pbmc
    count <- count + 1
    print(count)
    
  }
  names(each_sra) <- shrt
  #return(each_sra)
  if(length(SRA_in) > 1) {
  
    
    #
    object.features <- SelectIntegrationFeatures(object.list = each_sra, nfeatures = 3000)
    object.list <- PrepSCTIntegration(object.list = each_sra, anchor.features = object.features, 
                                        verbose = FALSE)
    #return(each_sra)
    # combine dataset with different integration anchors
    immune.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = object.features, normalization.method = "SCT")
    immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20, normalization.method = "SCT")
    
    DefaultAssay(immune.combined) <- "integrated"
    pbmc <- immune.combined
  }
  
  pbmc <- RunPCA(object = pbmc, verbose = FALSE)
  if(haveUmap == TRUE) {
    pbmc <- RunUMAP(object = pbmc, dims = 1:20, verbose = FALSE)
  }
  pbmc <- FindNeighbors(object = pbmc, dims = 1:20, verbose = FALSE)
  pbmc <- FindClusters(object = pbmc, verbose = FALSE)
  pbmc <- ScaleData(object = pbmc)
  
  if(saveALL == TRUE) {
    save(pbmc, file = paste0(name, "_custom.Rdata"))
  }
  return(pbmc)
}

process_dgTMatrix_lists <- function(dgTMatrix_list, name, species, naming_preference, panglao_set = FALSE ,haveUMAP = FALSE, saveSCObject = FALSE) {
sm <- dgTMatrix_list[[1]]
if(species == -9) {
  spec=get_gene_symbol(sm)
  species_name <- spec$species
}
pbmc <- process_from_count(countmat_list = dgTMatrix_list, name = name, theSpecies = species, panglao_set = panglao_set, haveUmap = haveUMAP, saveALL  = saveSCObject)
#return(pbmc)
print("A")
print(class(pbmc))
print(head(pbmc))
gsva_cellIdentity_out <- gsva_cellIdentify(pbmc, species_name, naming_preference)
save(gsva_cellIdentity_out, file = paste0(name, "gsva_cellname.Rdata"))
generes <- seurat_to_generes(pbmc = pbmc)
gene_out <- generes_to_heatmap(generes, name = name, species = species_name, naming_preference = naming_preference)

wilcoxon_map_rnk_t <- gene_out$pVal
wilcoxon_map_rnk_or <- gene_out$OR
names(generes) <- colnames(wilcoxon_map_rnk_t) <- colnames(wilcoxon_map_rnk_or) <- gene_out$cellname
save(generes, file = paste0(name, "_generes.Rdata"))
save(wilcoxon_map_rnk_t, file= paste0(name, "_pval_heatmap.Rdata"))
save(wilcoxon_map_rnk_or, file= paste0(name, "_or_heatmap.Rdata"))
return(wilcoxon_map_rnk_t)
}

get_naming_preference_options <- function() {
  load("cell_type_info/cell_preferences_categorized.Rdata")
  return(sort(names(cell_preference_final)))
}

get_gene_symbol <- function(wilcoxon_map_rnk_t) {
  the_human <- length(grep("ENSG00", rownames(wilcoxon_map_rnk_t)))
  the_mouse <- length(grep("ENSMUSG00", rownames(wilcoxon_map_rnk_t)))

  if(the_human >= the_mouse) {
    theSpecies <- "human"
  }
  if(the_mouse > the_human) {
    theSpecies <- "mouse"
  }
  sm <- wilcoxon_map_rnk_t
  
  # make sure that mitochondrial genes are flagged
  RN = rownames(sm)
  RN_1 <- sub('(.*)[.](.*)','\\1',RN)
  if(theSpecies == "human") {
    
    RN_2 <- substr(RN_1, 1, nchar(RN_1) - 16)
  }
  if(theSpecies == "mouse" ) {
    
    RN_2 <- substr(RN_1, 1, nchar(RN_1) - 19)
  }
  
  return(list(rowname = RN_2, species = theSpecies))
}

gsva_cellIdentify <- function(pbmc, theSpecies, naming_preference) {
  library(GSVA)
  avg_expr <- AverageExpression(pbmc)
  
  # panglao
  if(theSpecies == "human") {
    load("cell_type_info/human_marker_both.RData")
  } else {
    load("cell_type_info/mouse_marker_panglao.RData")
  }
  gbm <- gsva(as.matrix(avg_expr$RNA), gmt, mx.diff=FALSE, verbose=FALSE, parallel.sz=1, min.sz= 5)
  gbm_pang <- gbm[grep("panglao", rownames(gbm)),]
  gbm_cellmarker <- gbm[-grep("panglao", rownames(gbm)),]
  
  top5_pang <- list()
  top5_cm <- list()
  for(i in 1:ncol(gbm)) {
    top5_pang[[i]] <-  round(sort(gbm_pang[,i], decreasing = T)[1:5],2)
    top5_cm[[i]] <- round(sort(gbm_cellmarker[,i], decreasing = T)[1:5],2)
  }
  
  top_ct <- function(x) {
    if(naming_preference != -9) {
      load("cell_type_info/cell_preferences_categorized.Rdata")
      # if they have the prefered tissue or cell-type to pick from
      # prioritize those
      mypref <- cell_preference_final[[naming_preference]]
      top5_ct_pref <- names(x)
      outcell_preference <- which(top5_ct_pref %in% mypref)
      if(length(outcell_preference) > 0 ) {
        y <- x[outcell_preference[1]]
        return(paste0(names(y[1]),"_",unname(y[1])))
      } else {
        return(paste0(names(x[1]),"_",unname(x[1])))
      }
      
    } else {
      return(paste0(names(x[1]),"_",unname(x[1])))
    }
  }
  cm_top <- lapply(top5_cm, top_ct)
  pang_top <- lapply(top5_pang, top_ct)
  
  l <- list(cellMarker = cm_top, panglao = pang_top)
  return(l)
}
