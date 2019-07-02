##============================ Some supporting functions for functions for users ===========================##

## Check pheno variable
checkPheno <-  function(VAR, nameVAR)
{
  # Check whether VAR is a vector
  if(!is.vector(VAR)) {                             # If VAR is not a vector
    if(is.matrix(VAR) && any(dim(VAR)==1)) {        # Check whether VAR is a matrix of one row or one column
      VAR <- as.vector(VAR)                         # Convert the matrix into a vector
      warning( paste(nameVAR, "is a matrix!"), call. = FALSE)   # Also give warning
    }else
      stop(paste(nameVAR, "must be a vector."), call. = FALSE)  # Stop if VAR is neither a vector not a matrix of one row or column
  }
  # Check whether VAR is a numeric vector
  if(!is.numeric(VAR))
    stop(paste(nameVAR, "must be a numeric vector."), call. = FALSE)
  # Check whether there is any NA
  if(any(is.na(VAR)))
    stop(paste0("One or more phenotype values missing in ", nameVAR, "!"), call. = FALSE)
  # Check whether there is more than one non-missing arguments
  if(length(VAR) <= 1)
    stop(paste("Number of elements in the", nameVAR, "vector must be more than 1!"), call. = FALSE)
  return(VAR)
}

check_normalized_geno <- function(genoMbr, Index) {
  if(nrow(na.omit(genoMbr)) != nrow(genoMbr))           # If there is any missing observation
    stop(paste0("tissue", Index, " : Missing values in tissue-specific eQTL genotype matrix!"), call. = FALSE);
  if(all(as.vector(genoMbr) %% 1 == 0)) {               # If all the entries are integers
    if(all(unique(as.vector(genoMbr)) %in% 0:2)) {      # If all the entries are 0, 1 or 2
      frq <- colMeans(genoMbr)/2                        # Calculate allele frequency of all the SNPs
      if(any(frq > 0.99) || any(frq < 0.01))            # If MAF is very small
        stop(paste0("tissue", Index, " : MAF of each eQTL must be more than 1%."), call. = FALSE);
      genoMbr <- scale(genoMbr, center = TRUE, scale = TRUE)                                      # Otherwise normalize the genotype matrix
    }
    else stop(paste0("tissue", Index, " : eQTL genotypes must be coded as #minor allele; all eQTLs (SNPs) must be biallelic."), call. = FALSE)
  }else{
    ## check if the genotypes of an eQTL SNP are standardized/normalized.
    epsilon_dev_zero = 10^(-3)
    numCol = ncol(genoMbr)
    column_sds_dev = abs(matrixStats::colSds(genoMbr) - 1)
    column_means_dev = abs(colMeans(genoMbr))
    check_sds = sum(column_sds_dev < epsilon_dev_zero)     ## check how many sds = 1
    check_means = sum(column_means_dev < epsilon_dev_zero)     ## check how many means = 0
    ## any mean non-zero ?
    if(check_means < numCol)
      stop(paste0("tissue", Index, " : eQTL genotypes not normalized, non-zero mean"), call. = FALSE);
    ## any sd non-unity ?
    if(check_sds < numCol)
      stop(paste0("tissue", Index, " : eQTL genotypes not normalized, non-unity variance"), call. = FALSE);
  }
  return(genoMbr)
}

check_geno <- function(geno, pheno) {
  if(length(geno) == 1 && is.na(geno))                    # If geno is missing
    stop("geno is missing!", call. = FALSE)
  if(!is.list(geno))                                      # If geno is not a list
    stop("geno must be a list.", call. = FALSE)
  if(length(geno) < 2)                                    # If geno is a list, but there are less than 2 elements
    stop("geno must be a list containing at least two matrices.", call. = FALSE)
  for(i in 1:length(geno)) {                              # Now check each element of the list
    if(length(geno[[i]]) == 1 && is.na(geno[[i]]))
      stop("Element of geno is missing", call. = FALSE)
    if(is.data.frame(geno[[i]])){                         # If ith element is a dataframe
      geno[[i]] <- as.matrix(geno[[i]])                   # Convert it into a matrix
      if(!is.numeric(geno[[i]]))                          # If the converted matrix is not numeric
        stop("Element of geno is not numeric!", call. = FALSE)  # Stop
      warning("Element of geno is data.frame, being converted to matrix.", call. = FALSE) # If converted matrix is numeric, let user know about the conversion
    }
    else if(is.vector(geno[[i]]))
      stop("Element of geno is not a matrix. It's vector.", call. = FALSE)
    else if(!is.matrix(geno[[i]]) || !is.numeric(geno[[i]]))   # If ith element is either not a matrix or not numeric
      stop("Element of geno is not a numeric matrix!", call. = FALSE) # Stop
    if(nrow(geno[[i]]) != length(pheno))                  # If number of rows of ith element matrix is not equal to number of subjects
      stop("Number of rows of genotype matrix must be the same as length of phenotype vector!", call. = FALSE)
    geno[[i]] <- check_normalized_geno(geno[[i]], i)      # Check whether the matrix is normalized
  }
  return(geno)
}

## Check input variable 'tissues'
check_tissues <-  function(tissues, geno)
{
  # Check whether argument 3 is a vector
  if(!is.vector(tissues))
    stop("tissues must be a vector!", call. = FALSE)
  # Check whether argument 3 is a character vector
  if(!is.character(tissues))
    stop("tissues must be a character vector!", call. = FALSE)
  # Check whether there is duplicate tissues
  if(length(tissues) > length(unique(tissues)))
    stop("Two or more tissues have the same name!", call. = FALSE)
  # Check whether argument 3 is a vector of length more than 1
  if(length(tissues) != length(geno))
    stop("Number of tissues and elements in list geno must be the same!", call. = FALSE)
}


