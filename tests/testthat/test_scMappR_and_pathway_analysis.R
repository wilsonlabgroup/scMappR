data(PBMC_example)
bulk_DE_cors <- PBMC_example$bulk_DE_cors
bulk_normalized <- PBMC_example$bulk_normalized
odds_ratio_in <- PBMC_example$odds_ratio_in
max_proportion_change <- 10
print_plots <- FALSE

testthat::test_that('Bulk samples, signature matrices, correct DEGs and index names.', {
  
case_grep <- "_female"
control_grep <- "_male"

toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
                                      bulk_DE_cors, case_grep = case_grep,
                                      control_grep = control_grep, rda_path = "", 
                                      max_proportion_change = 10, print_plots = TRUE, 
                                      plot_names = "tst1", theSpecies = "human", 
                                      output_directory = "tester",
                                      sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                      internet = FALSE)
})


testthat::test_that('Must have multiple replicates -- throw error', {
  
  case_grep <- "_female"
  control_grep <- "_male"
  bulk_normalized1 <- bulk_normalized[,1:4]
  
  
  toOut <- scMappR_and_pathway_analysis(bulk_normalized1, odds_ratio_in, 
                                        bulk_DE_cors, case_grep = case_grep,
                                        control_grep = control_grep, rda_path = "", 
                                        max_proportion_change = 10, print_plots = TRUE, 
                                        plot_names = "tst1", theSpecies = "human", 
                                        output_directory = "tester",
                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                        internet = FALSE)
})


testthat::test_that('Must have multiple replicates -- throw error', {
  
  case_grep <- "_FEMALE"
  control_grep <- "_MALE"

  
  
  toOut <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
                                        bulk_DE_cors, case_grep = case_grep,
                                        control_grep = control_grep, rda_path = "", 
                                        max_proportion_change = 10, print_plots = TRUE, 
                                        plot_names = "tst1", theSpecies = "human", 
                                        output_directory = "tester",
                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                        internet = FALSE)
})