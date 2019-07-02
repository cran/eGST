context("Input EQTL genotype data")

test_that("Check whether EQTL genotype data is a list of two matrices", {
  expect_error(eGST(ExamplePhenoData), "geno is missing!")
  expect_error(eGST(ExamplePhenoData, geno = NA, paste0("tissue", 1:2)), "geno is missing!")
  expect_error(eGST(ExamplePhenoData, geno = 1:2, paste0("tissue", 1:2)), "geno must be a list.")
  expect_error(eGST(ExamplePhenoData, geno = matrix(1:10, 2, 5), paste0("tissue", 1:2)), "geno must be a list.")
  expect_error(eGST(ExamplePhenoData, geno = list(matrix(1:10, 2, 5)), paste0("tissue", 1:2)), "geno must be a list containing at least two matrices.")
  expect_error(eGST(ExamplePhenoData, geno = list(data.frame(1:10, 2, 5)), paste0("tissue", 1:2)), "geno must be a list containing at least two matrices.")
  expect_error(eGST(ExamplePhenoData, geno = list(ExampleEQTLgenoData[[1]][,1], ExampleEQTLgenoData[[2]][,1]), paste0("tissue", 1:2)), "Element of geno is not a matrix. It's vector.")
  expect_error(eGST(ExamplePhenoData, geno = data.frame(1:10, 2, 5), paste0("tissue", 1:2)), "Element of geno is not a matrix. It's vector.")
  expect_error(eGST(ExamplePhenoData, geno = list(1:2, 3:4), paste0("tissue", 1:2)), "Element of geno is not a matrix. It's vector.")
  expect_warning(eGST(ExamplePhenoData, geno = list(as.data.frame(ExampleEQTLgenoData[[1]]), as.data.frame(ExampleEQTLgenoData[[2]])), paste0("tissue", 1:2)), "Element of geno is data.frame, being converted to matrix.")
})

test_that("Check whether elements of geno are numeric matrix", {
  expect_error(eGST(ExamplePhenoData, geno = list(ExampleEQTLgenoData[[1]], ExampleEQTLgenoData[[2]][1,]), paste0("tissue", 1:2)), "Element of geno is not a matrix. It's vector.")
  expect_error(eGST(ExamplePhenoData, geno = list(matrix(letters[1:10], 5, 2), matrix(letters[1:10], 5, 2)), paste0("tissue", 1:2)), "Element of geno is not a numeric matrix!")
  expect_error(eGST(ExamplePhenoData, geno = list(as.data.frame(matrix(letters[1:10], 5, 2)), ExampleEQTLgenoData[[1]]), paste0("tissue", 1:2)), "Element of geno is not numeric!")
})

test_that("Check whether each genotype matrix has same number of rows as the length of phenotype vector", {
  expect_error(eGST(ExamplePhenoData, geno = list(ExampleEQTLgenoData[[1]], ExampleEQTLgenoData[[2]][1:10,]), paste0("tissue", 1:2)), "Number of rows of genotype matrix must be the same as length of phenotype vector!")
})

test_that("Check whether there is any missing entry in genotype matrix", {
  expect_error(eGST(ExamplePhenoData, list(matrix(sample(c(0:2, NA), length(ExamplePhenoData)*10, replace = TRUE), length(ExamplePhenoData)), ExampleEQTLgenoData[[1]]), paste0("tissue", 1:2)), "tissue1 : Missing values in tissue-specific eQTL genotype matrix!")
  expect_error(eGST(ExamplePhenoData, list(matrix(sample(c(0:2), length(ExamplePhenoData)*10, replace = TRUE), length(ExamplePhenoData)), NA), paste0("tissue", 1:2)), "Element of geno is missing")
})

test_that("Check non-normalized genotype matrix", {
  expect_error(eGST(ExamplePhenoData, list(matrix(sample(c(1:3), length(ExamplePhenoData)*10, replace = TRUE), length(ExamplePhenoData)), ExampleEQTLgenoData[[1]]), paste0("tissue", 1:2)), "tissue1 : eQTL genotypes must be coded as *")
  expect_error(eGST(ExamplePhenoData, list(matrix(sample(c(0.1:0.3), length(ExamplePhenoData)*10, replace = TRUE), length(ExamplePhenoData)), ExampleEQTLgenoData[[1]]), paste0("tissue", 1:2)), "tissue1 : eQTL genotypes not normalized, non-zero mean")
  expect_error(eGST(ExamplePhenoData, list(cbind(rep(0, length(ExamplePhenoData)), sample(0:2, length(ExamplePhenoData), replace = TRUE)), ExampleEQTLgenoData[[1]]), paste0("tissue", 1:2)), "tissue1 : MAF of each eQTL must be more than 1%.")
})




