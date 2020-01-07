

testthat::test_that('60 genes within the matrix for POA. Correct dim and class', {
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  sig <- get_gene_symbol(POA_Rank_signature)
  Signature <- POA_Rank_signature
  rownames(Signature) <- sig$rowname
genes <- rownames(Signature)[1:60]
heatmap_test <- tissue_scMappR_custom( genes, signature_matrix = Signature,
                                       output_directory =  "scMappR_test", toSave = FALSE)
testthat::expect_equal(dim(heatmap_test$gene_list_heatmap$geneHeat), c(60,27))
testthat::expect_true(is.list(heatmap_test))
})

testthat::test_that('One gene, Error due to downstream analysis.', {
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  sig <- get_gene_symbol(POA_Rank_signature)
  Signature <- POA_Rank_signature
  rownames(Signature) <- sig$rowname
  genes <- rownames(Signature)[1:1]
  testthat::expect_error(tissue_scMappR_custom( genes, signature_matrix = Signature,
                                         output_directory =  "scMappR_test", toSave = FALSE))
})

testthat::test_that('Three genes, Small list should work.', {
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  sig <- get_gene_symbol(POA_Rank_signature)
  Signature <- POA_Rank_signature
  rownames(Signature) <- sig$rowname
  genes <- rownames(Signature)[1:3]
  heatmap_test <- tissue_scMappR_custom( genes, signature_matrix = Signature,
                                         output_directory =  "scMappR_test", toSave = FALSE)
  
  testthat::expect_equal(dim(heatmap_test$gene_list_heatmap$geneHeat), c(3,27))
  testthat::expect_true(is.list(heatmap_test))
})

testthat::test_that('Mix of genes in signature and other characters', {
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  sig <- get_gene_symbol(POA_Rank_signature)
  Signature <- POA_Rank_signature
  rownames(Signature) <- sig$rowname
  genes <- c(rownames(Signature)[1:3], 1, "2", "a")
  heatmap_test <- tissue_scMappR_custom( genes, signature_matrix = Signature,
                                         output_directory =  "scMappR_test", toSave = FALSE)
  
  testthat::expect_equal(dim(heatmap_test$gene_list_heatmap$geneHeat), c(3,27))
  testthat::expect_true(is.list(heatmap_test))
})


