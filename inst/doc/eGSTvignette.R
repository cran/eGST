## ----install_package, eval=FALSE, collapse = TRUE------------------------
#  #install.packages("eGST")
#  #library("eGST")

## ----load_pheno, collapse = TRUE-----------------------------------------
library("eGST")
# Load the phenotype data vector
phenofile <- system.file("extdata", "ExamplePhenoData.rda", package = "eGST")
load(phenofile)
head(ExamplePhenoData)

## ----load_geno, collapse = TRUE------------------------------------------
library("eGST")
# Load the list containing genotype matrices of tissue-specific eQTLs. 
genofile <- system.file("extdata", "ExampleEQTLgenoData.rda", package = "eGST")
load(genofile)
ExampleEQTLgenoData[[1]][1:5, 1:5]

## ----names, collapse = TRUE----------------------------------------------
# Specify the name of the tissues.
tissues <- paste0("tissue", 1:2)
tissues

## ----example_run_eGST, collapse = TRUE-----------------------------------
#Run eGST to estimate the tissue-specific posterior probability across individuals.
result <- eGST(ExamplePhenoData, ExampleEQTLgenoData, tissues, nIter = 10)

## ----result_structure, collapse= TRUE------------------------------------
# Overall summary of the results produced by eGST.
str(result)

