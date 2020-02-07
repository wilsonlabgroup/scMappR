testthat::test_that('Toy dataset, should run', {
bulk_DE_cors <- cbind(1:50, stats::runif(n = 50, min = 0, 0.05), stats::runif(50, -5, 5))
colnames(bulk_DE_cors) <- c("gene_name", "padj", "log2fc")
rownames(bulk_DE_cors) <- 1:50
bulk_normalized <- matrix(stats::runif(2000*10,0, 100), ncol = 10)
colnames(bulk_normalized) <- c(paste0(1:5,"_female"), paste0(1:5,"_male"))
rownames(bulk_normalized) <- 1:2000
odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
rownames(odds_ratio_in) <- c(30:129)
colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
max_proportion_change <- 10
print_plots <- FALSE  
case_grep <- "_female"
control_grep <- "_male"

tst1 <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
                             bulk_DE_cors, case_grep = case_grep,
                             control_grep = control_grep, rda_path = "", 
                             max_proportion_change = 10, print_plots = TRUE, 
                             plot_names = "tst1", theSpecies = "human", 
                             output_directory = "tester",
                             sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                             internet = FALSE)
testthat::expect_true(all(names(tst) == c("cellWeighted_Foldchange","cellType_Proportions","leave_one_out_proportions","processed_signature_matrix","ProportionT.test")))
})

testthat::test_that('Must have multiple replicates -- throw error', {
  bulk_DE_cors <- cbind(1:50, stats::runif(n = 50, min = 0, 0.05), stats::runif(50, -5, 5))
  colnames(bulk_DE_cors) <- c("gene_name", "padj", "log2fc")
  rownames(bulk_DE_cors) <- 1:50
  bulk_normalized <- matrix(stats::runif(2000*10,0, 100), ncol = 10)
  colnames(bulk_normalized) <- c(paste0(1:5,"_female"), paste0(1:5,"_male"))
  rownames(bulk_normalized) <- 1:2000
  odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
  rownames(odds_ratio_in) <- c(30:129)
  colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
  max_proportion_change <- 10
  print_plots <- FALSE  
  case_grep <- "_female"
  control_grep <- "_male"
  bulk_normalized1 <- bulk_normalized[,1:4]
  
  
  toOut <- testthat::expect_error(scMappR_and_pathway_analysis(bulk_normalized1, odds_ratio_in, 
                                        bulk_DE_cors, case_grep = case_grep,
                                        control_grep = control_grep, rda_path = "", 
                                        max_proportion_change = 10, print_plots = TRUE, 
                                        plot_names = "tst1", theSpecies = "human", 
                                        output_directory = "tester",
                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                        internet = FALSE))
})


testthat::test_that('Strange column names, same error', {
  
  bulk_DE_cors <- cbind(1:50, stats::runif(n = 50, min = 0, 0.05), stats::runif(50, -5, 5))
  colnames(bulk_DE_cors) <- c("gene_name", "padj", "log2fc")
  rownames(bulk_DE_cors) <- 1:50
  bulk_normalized <- matrix(stats::runif(2000*10,0, 100), ncol = 10)
  colnames(bulk_normalized) <- c(paste0(1:5,"_female"), paste0(1:5,"_male"))
  rownames(bulk_normalized) <- 1:2000
  odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
  rownames(odds_ratio_in) <- c(30:129)
  colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
  max_proportion_change <- 10
  print_plots <- FALSE 
  case_grep <- "_FEMALE"
  control_grep <- "_MALE"

  
  toOut <- testthat::expect_error(scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
                                        bulk_DE_cors, case_grep = case_grep,
                                        control_grep = control_grep, rda_path = "", 
                                        max_proportion_change = 10, print_plots = TRUE, 
                                        plot_names = "tst1", theSpecies = "human", 
                                        output_directory = "tester",
                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                        internet = FALSE))
})

testthat::test_that('Case and control indexed the same, throw error.', {
  
  bulk_DE_cors <- cbind(1:50, stats::runif(n = 50, min = 0, 0.05), stats::runif(50, -5, 5))
  colnames(bulk_DE_cors) <- c("gene_name", "padj", "log2fc")
  rownames(bulk_DE_cors) <- 1:50
  bulk_normalized <- matrix(stats::runif(2000*10,0, 100), ncol = 10)
  colnames(bulk_normalized) <- c(paste0(1:5,"_female"), paste0(1:5,"_male"))
  rownames(bulk_normalized) <- 1:2000
  odds_ratio_in <- matrix(stats::runif(100*5,0,6),ncol=5 )
  rownames(odds_ratio_in) <- c(30:129)
  colnames(odds_ratio_in) <- paste0("celltype_", 1:5)
  max_proportion_change <- 10
  print_plots <- FALSE  
  case_grep <- c(1,2,3,4)
  control_grep <- "_female"
  
  
  
  toOut <- testthat::expect_error(scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in, 
                                        bulk_DE_cors, case_grep = case_grep,
                                        control_grep = control_grep, rda_path = "", 
                                        max_proportion_change = 10, print_plots = TRUE, 
                                        plot_names = "tst1", theSpecies = "human", 
                                        output_directory = "tester",
                                        sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                        internet = FALSE))
})

testthat::test_that('one cell type, throw error.', {
  
  
  data(PBMC_example)
  bulk_DE_cors <- PBMC_example$bulk_DE_cors
  bulk_normalized <- PBMC_example$bulk_normalized
  odds_ratio_in <- PBMC_example$odds_ratio_in
  max_proportion_change <- 10
  print_plots <- FALSE  
  case_grep <- "_female"
  control_grep <- "_male"
  
  
  
  toOut <- testthat::expect_error(scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in[,1], 
                                                               bulk_DE_cors, case_grep = case_grep,
                                                               control_grep = control_grep, rda_path = "", 
                                                               max_proportion_change = 10, print_plots = TRUE, 
                                                               plot_names = "tst1", theSpecies = "human", 
                                                               output_directory = "tester",
                                                               sig_matrix_size = 3000, up_and_downregulated = FALSE, 
                                                               internet = FALSE))
})