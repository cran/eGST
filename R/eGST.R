#' eGST (eQTL-based Genetic Sub-Typer): Leveraging eQTLs to identify
#'  individual-level tissue of interest for a complex trait.
#'
#'
#'Genetic predisposition for complex traits is often manifested through multiple tissues
#' of interest at different time points in the development. As an example, the genetic predisposition
#'  for obesity could be manifested through inherited variants that control metabolism through regulation
#'   of genes expressed in the brain and/or through the control of fat storage in the adipose tissue by
#'    dysregulation of genes expressed in adipose tissue. We present a method eGST that integrates
#'     tissue-specific eQTLs with GWAS data for a complex trait to probabilistically assign a
#'      tissue of interest to the phenotype of each individual in the study. eGST estimates
#'       the posterior probability that an individual's phenotype can be assigned to a
#'        tissue based on individual-level genotype data of tissue-specific eQTLs and marginal
#'         phenotype data in a GWAS cohort. Under a Bayesian framework of mixture model,
#'          eGST employs a maximum a posteriori (MAP) expectation-maximization (EM) algorithm
#'           to estimate the tissue-specific posterior probability across individuals.
#' @aliases eGST-package
#'
# \code{\link{eGST-package}}
#'
#' @section Functions:
#' \describe{
#' \item{\code{\link{eGST}}}{It estimates the posterior probability that the genetic susceptibility
#'  of the phenotype of an individual in the study is mediated through eQTLs specific to a tissue
#'   of interest. The phenotype across individuals can be classified into tissues under consideration
#'    based on the estimated tissue-specific posterior probability across individuals.}}
#' @references Majumdar A, Giambartolomei C, Cai N, Freund MK, Haldar T, J Flint, Pasaniuc B (2019) Leveraging eQTLs to identify
#' tissue-specific genetic subtype of complex trait. bioRxiv.
#'
#' @docType package
#'
#'
#' @importFrom matrixStats colSds
#' @importFrom mvtnorm dmvnorm
#' @importFrom purrr map_dbl map_int
#' @importFrom MASS mvrnorm
#' @importFrom utils combn
#' @importFrom stats runif rnorm rbeta quantile qchisq qbeta pchisq pbeta p.adjust dist aggregate sd dnorm qnorm var na.omit
"_PACKAGE"
