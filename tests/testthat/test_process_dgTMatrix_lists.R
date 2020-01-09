
testthat::test_that('cell-type to convert to signature matrix should work', {
  
data(sm)
toProcess <- list(example = sm)
tst1 <- process_dgTMatrix_lists(toProcess, name = "testPropcess", species_name = -9,
                                naming_preference = "eye", rda_path = "~/scMappR/data", panglao_set = TRUE)

output_list <- all(names(tst1) == c("wilcoxon_rank_mat_t", "wilcoxon_rank_mat_or", "generes"))
testthat::expect_true(output_list)
testthat::expect_true(is.list(tst1))

})

testthat::test_that('too few cells, should fail', {
  
  data(sm)
  toProcess <- list(example = sm[,1:20])
   testthat::expect_error(process_dgTMatrix_lists(toProcess, name = "testPropcess", species_name = -9,
                                  naming_preference = "eye", rda_path = "~/scMappR/data", panglao_set = TRUE))
})
