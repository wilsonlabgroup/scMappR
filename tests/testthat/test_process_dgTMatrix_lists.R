

testthat::test_that('too few cells, should fail', {
  
  data(sm)
  toProcess <- list(example = sm[,1:20])
   testthat::expect_error(process_dgTMatrix_lists(toProcess, name = "testPropcess", species_name = -9,
                                  naming_preference = "eye", rda_path = "~/scMappR/data"))
})
