testthat::test_that('Twenty genes should work', {
  odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
  rownames(odds_ratio_in) <- c(30:129)
  colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
  Signature <- odds_ratio_in
  genes <- rownames(Signature)[1:20]
  tst1 <- tissue_scMappR_custom( genes, signature_matrix = Signature,
                                                output_directory =  "scMappR_test", toSave = FALSE)
  testthat::expect_true(all(colnames(tst1) == c("background_heatmap", "gene_list_heatmap", "single_celltype_preferences","group_celltype_preferences")  ))
  
})




testthat::test_that('One gene, Error due to downstream analysis.', {
  odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
  rownames(odds_ratio_in) <- c(30:129)
  colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
  Signature <- odds_ratio_in
  
  genes <- rownames(Signature)[1]
  testthat::expect_error(tissue_scMappR_custom( genes, signature_matrix = Signature,
                                         output_directory =  "scMappR_test", toSave = FALSE))
})



