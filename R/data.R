# Documentation of the example dataset to demonstrate how to run eGST.
#'An example of phenotype data.
#'
#'ExamplePhenoData is a simulated vector of phenotype for 1000 individuals.
#'In this simulated example dataset, we have considered two tissues and corresponding
#'sets of 100 tissue-specific eQTLs each.
#'First half of 1000 individuals' phenotypes were simulated to have genetic effect from the first tissue
#'specific eQTLs, but no effect from the second tissue-specific eQTLs. Hence first 500 individuals
#'were assigned to the first tissue. Similarly, second half of
#'the 1000 individuals were simulated to have genetic effect from the second-tissue specific eQTLs and hence
#'they were assigned to the second tissue.
#'
#'@usage data(ExamplePhenoData)
#'
#'@format A numeric vector of length 1000 containing the phenotype data of 1000 individuals.
#'
#'@examples
#'data(ExamplePhenoData)
#'pheno <- ExamplePhenoData
"ExamplePhenoData"

#'An example of tissue-specific eQTLs genotype data.
#'
#'ExampleEQTLgenoData is a list with two elements each containing a 1000 by 100 ordered genotype matrix.
#'Each matrix provides the genotype data of the 1000 individuals at 100 tissue-specific eQTLs for each tissue.
#'
#'@usage data(ExampleEQTLgenoData)
#'
#'@format A list of two numeric matrix each having 1000 rows (individuals) and 100 columns (eQTLs):
#'
#'@examples
#'data(ExampleEQTLgenoData)
#'geno <- ExampleEQTLgenoData
"ExampleEQTLgenoData"
