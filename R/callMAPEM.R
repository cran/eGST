## Function to call the main MAP-EM function for implementing eGST.
#' Run eGST.
#'
#' Run eGST to estimate the posterior probability that the genetic
#' susceptibility of the phenotype of
#' an individual in the study is mediated through eQTLs specific to a tissue of interest.
#' To create sets of tissue-specific eQTLs in your context, please see our manuscript:
#'  Majumdar A, Giambartolomei C, Cai N, Freund MK, Haldar T, J Flint, Pasaniuc B (2019)
#'   Leveraging eQTLs to identify tissue-specific genetic subtype of complex trait, bioRxiv.
#' @param pheno A numeric vector of length N where N is the number of individuals. It
#'  contains the GWAS phenotype values of individuals. No default.
#' @param geno A list with K elements where K is the number of tissues. Each element of
#' geno is the genotype matrix of the eQTLs specific to a tissue in the GWAS cohort.
#'  So j-th element of geno is N by Mj matrix
#' containing the genotype data of N individuals (rows) at Mj eQTLs (columns) specific to j-th tissue.
#' Each eQTL is a bi-allelic SNP with minor allele frequency > 0.01.
#' Genotypes at each eQTL must be normalized across N individuals.
#' If 0/1/2 valued genotype matrix is provided, it is internally normalized. No default.
#' @param tissues A character vector of length K. It contains the names of tissues of interest in the analysis.
#'  The order of tissues in this vector must match the order of tissues in the previous argument 'geno'. No default.
#' @param logLimprovement A positive real number specifying the minimum possible improvement
#'  of data log-likelihood in MAP-EM stopping criterion. Default \eqn{5*10^(-8)}.
#' @param seed_choice An integer providing the choice of random seed for initialization in MAP-EM algorithm.
#'  Default is an integer randomly selected in (1,...,1000).
#' @param nIter An integer providing the maximum number of iterations allowed in the MAP-EM algorithm. Default is 100.
#' @return The output produced by \code{\link{eGST}} is a list which consists of various components.
#'    \item{gamma}{A N by K matrix providing the tissue-specific posterior probability of N individuals across K tissues.}
#'    \item{alfa}{Baseline tissue-specific intercepts/means of the trait.}
#'    \item{beta}{Tissue-specific eQTLs' genetic effect on the trait.}
#'    \item{sigma_g}{Square root of the variance of tissue-specific per-eQTL genetic effect on the trait.}
#'    \item{sigma_e}{Square root of the error variance of tissue-specific subtype of the trait
#'     which remains unexplained by the tissue-specific eQTLs.}
#'    \item{m}{Number of tissue-specific eQTLs.}
#'    \item{logL}{log-likelihood of the data.}
#'
#' @references A Majumdar, C Giambartolomei, N Cai, MK Freund, T Haldar, T Schwarz, J Flint, B Pasaniuc (2019) Leveraging eQTLs to identify
#' tissue-specific genetic subtype of complex trait. bioRxiv.
#'
#'
#' @examples
#' data(ExamplePhenoData)
#' pheno <- ExamplePhenoData
#' head(pheno)
#' data(ExampleEQTLgenoData)
#' geno <- ExampleEQTLgenoData
#' geno[[1]][1:5,1:5]
#' geno[[2]][1:5,1:5]
#' tissues <- paste("tissue", 1:2, sep = "")
#' result <- eGST(pheno, geno, tissues)
#' str(result)
#'
#' @export
eGST <- function(pheno, geno, tissues, logLimprovement = 5*10^(-8), seed_choice = sample(1:1000, size=1), nIter = 100)
{
  ## check if any of the main arguments is missing.
  if(missing(pheno)) stop("pheno vector is missing!", call. = FALSE)
  if(missing(geno)) stop("geno is missing!", call. = FALSE)
  if(missing(tissues)) stop("tissues vector is missing!", call. = FALSE)
  ## check the pheno vector, a numeric vector
  pheno <- checkPheno(pheno, "pheno")
  ## Check geno list
  geno <- check_geno(geno, pheno)
  ## check if tissues and geno have desirable properties.
  check_tissues(tissues, geno)
  nTissues = length(tissues);

    # Rename the arguments and call the MAPEM function.
    Y = pheno
    X = geno
    REPLI = nIter

    #invisible(capture.output(RESULT <- MAPEM(Y, X, tissues, nTissues, logLimprovement, seed_choice, REPLI)))
    RESULT <- MAPEM(Y, X, tissues, nTissues, logLimprovement, seed_choice, REPLI)
    #print_result(RESULT)
    #invisible(RESULT)
}



