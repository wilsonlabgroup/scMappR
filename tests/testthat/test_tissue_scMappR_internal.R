
testthat::test_that('200 genes hypothal should work', {
  
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  Signature <- POA_Rank_signature
  rowname <- get_gene_symbol(Signature)
  rownames(Signature) <- rowname$rowname
  
  rda_path1 = "~/scMappR/data"
  
genes <- rownames(Signature)[1:200]
internal <- tissue_scMappR_internal(genes,"mouse", output_directory = "scMappR_Test",
                                    tissue = "hypothalamus",rda_path = rda_path1)

CNames <- c("background_heatmap","gene_list_heatmap" ,"single_celltype_preferences","group_celtype_preference"   )

allTrue <- all(names(internal[[1]]) == CNames)

testthat::expect_true(allTrue)
testthat::expect_true(is.list(internal))

})


testthat::test_that('unavailable cell-type should fail', {
  
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  Signature <- POA_Rank_signature
  rowname <- get_gene_symbol(Signature)
  rownames(Signature) <- rowname$rowname
  
  rda_path1 = "~/scMappR/data"
  
  genes <- rownames(Signature)[1:200]
  testthat::expect_error(tissue_scMappR_internal(genes,"mouse", output_directory = "scMappR_Test",
                                      tissue = "potato",rda_path = rda_path1))
  
})


testthat::test_that('wrong species should throw error', {
  
  data(POA_example)
  POA_generes <- POA_example$POA_generes
  POA_OR_signature <- POA_example$POA_OR_signature
  POA_Rank_signature <- POA_example$POA_Rank_signature
  Signature <- POA_Rank_signature
  rowname <- get_gene_symbol(Signature)
  rownames(Signature) <- rowname$rowname
  
  rda_path1 = "~/scMappR/data"
  
  genes <- rownames(Signature)[1:200]
  internal <- testthat::expect_error(tissue_scMappR_internal(genes,"human", output_directory = "scMappR_Test",
                                      tissue = "hypothalamus",rda_path = rda_path1))
  

  
})

