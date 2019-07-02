context("Input phenotype data")

test_that("Check whether phenotype data is a numeric vector", {
  expect_error(eGST(), "pheno vector is missing!")
  expect_error(eGST("1", ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno must be a numeric vector.")
  expect_error(eGST(matrix(1:10, 2, 5), ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno must be a vector.")
  expect_error(eGST(matrix(letters[1:10], 2, 5), ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno must be a vector.")
  expect_warning(eGST(matrix(ExamplePhenoData, 1, 1000), ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno is a matrix!")
  expect_warning(eGST(matrix(ExamplePhenoData, 1000, 1), ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno is a matrix!")
  expect_error(eGST(data.frame(a = 1:10, b = letters[1:10]), ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno must be a vector.")
  expect_error(eGST(list(1:10), ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno must be a numeric vector.")
})

test_that("Check for missing values in the phenotype data", {
  expect_error(eGST(c(1:10, NA), ExampleEQTLgenoData, paste0("tis", 1:2)), "One or more phenotype values missing in pheno!")
  expect_error(eGST(NA, ExampleEQTLgenoData, paste0("tis", 1:2)), "pheno must be a numeric vector.")
})

test_that("Check whether phenotype data contains more than one element", {
  expect_error(eGST(1, ExampleEQTLgenoData, paste0("tis", 1:2)), "Number of elements in the pheno vector must be more than 1!")
})
