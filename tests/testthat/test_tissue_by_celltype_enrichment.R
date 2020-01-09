

testthat::test_that('100 genes hypothal should work', {
  
  
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  Signature <- POA_Rank_signature
  rowname <- get_gene_symbol(Signature)
  rownames(Signature) <- rowname$rowname
  genes <- rownames(Signature)[1:100]
  
enriched <- tissue_by_celltype_enrichment(gene_list = genes, 
                                          species = "mouse",p_thresh = 0.05, isect_size = 3)

CName <- c("name", "p", "term_size", "intersect_size",
           "input_length", "genes", "fdr", "bonf", 
           "term_name","p_value")

testthat::expect_true(all(colnames(enriched) == CName))


})

testthat::test_that('Wrong species should throw warning and have no enriched cell-types', {
  
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  Signature <- POA_Rank_signature
  rowname <- get_gene_symbol(Signature)
  rownames(Signature) <- rowname$rowname
  genes <- rownames(Signature)[1:100]
  enriched <- tissue_by_celltype_enrichment(gene_list = genes, 
                                            species = "human",p_thresh = 0.05, isect_size = 3)
 testthat::expect_true(all(dim(enriched) == c(0, 10)))
  
})

