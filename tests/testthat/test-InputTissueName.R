context("Check input vector of tissue name")

test_that("Check whether tissue name is a character vector", {
  expect_error(eGST(ExamplePhenoData, ExampleEQTLgenoData, matrix(1:2, 1, 2)), "tissues must be a vector!")
  expect_error(eGST(ExamplePhenoData, ExampleEQTLgenoData, 1:2), "tissues must be a character vector!")
  expect_error(eGST(ExamplePhenoData, ExampleEQTLgenoData, list(1)), "tissues must be a character vector!")
  expect_error(eGST(ExamplePhenoData, ExampleEQTLgenoData, "1"), "Number of tissues and elements in list geno must be the same!")
  expect_error(eGST(ExamplePhenoData, ExampleEQTLgenoData, c("a", "a")), "Two or more tissues have the same name!")
})
