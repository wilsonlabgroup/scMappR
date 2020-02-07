


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



