## This function runs the 2 component EM for sub-phenotyping

MAPEM = function(Y, X, tissues, nTissues, logLimprovement, seed_choice, REPLI){

  #suppressMessages(library("mvtnorm"))
  #suppressMessages(library("MASS"))
  #suppressMessages(library("purrr"))

  #source("MAPEM_functions.R")               #### Source the core MAP-EM functions

  message("-------MAPEM in eGST starting--------")
  message(Sys.time())

  set.seed(seed_choice)

  m = purrr::map_int(X, ncol)       ## number of SNPs
  N = length(Y)              ## number of individuals

  ###================== Initialize error variances ====================##

  sigma_e = rep( sqrt( 0.95 * var(Y) ), nTissues )

  ###================== Initialize genetic variances ==================##

  sigma_g = numeric(nTissues)
  for(j in 1:nTissues) sigma_g[j] = sqrt( ( 0.05 * var(Y) ) / m[j] )

  ###================= Initialize alfa, beta ==================###

  alfa = rep(0, nTissues); sigma_alfa0 = sqrt(0.1);

  beta = as.list(1:nTissues)

  for(j in 1:nTissues) beta[[j]] = rnorm( m[j], mean = 0, sd = sigma_g[j] )

  ##----------------- Specify hyper parameters ----------------##

  a_e = 4                              ## sigE ~ Inverse-Gamma(a_e, b_e)
  b_e = (a_e - 1) * 0.95 * var(Y)      ## assuming 95% of total phenotype variance being explained by the error part

  a_g = 4
  b_g = (a_g - 1) * 0.05 * var(Y)      ## assuming 5% of total phenotype variance being explained by the genetic part

  ##----- Initialize the probability of the two sub-classes -------#

  w = rep(1/nTissues, nTissues)

  ## Initialize gamma: posterior probability of different phenotype observations to subclasses

  results = optim_gamma(N, nTissues, Y, X, alfa, beta, w, sigma_e)
  old_logL = results$logL

  ## starting the EM algorithm loop
  repli = 0;
  improvement = 10; logL_epsilon = logLimprovement;     ## related to log-likelihood change of the full data likelihood


  while(improvement > logL_epsilon && repli < REPLI){

    repli = repli+1
    message(paste0("Iteration ", repli, ":"))

    ##============== E-step ================##

    results = optim_gamma(N, nTissues, Y, X, alfa, beta, w, sigma_e)
    gamma = results$gamma
    new_logL = results$logL

    n = optim_n(gamma)

    ##=============== M-step ================##

    w = optim_w(n)

    ## update alfa
    alfa = optim_alfa( nTissues, N, n, gamma, Y, X, beta, sigma_e, sigma_alfa0 )

    ## update beta
    for(j in 1:nTissues) beta[[j]] = optim_beta( j, m, N, X, gamma, alfa, Y, sigma_e, sigma_g )

    ## update sigma_e
    sigma_e = optim_sigma_e( nTissues, N, Y, gamma, X, alfa, beta, n, a_e, b_e )

    ## update sigma_g
    sigma_g = optim_sigma_g( nTissues, m, beta, sigma_e, a_g, b_g )

    improvement = new_logL - old_logL
    if(repli == 1) improvement = 10
    old_logL = new_logL

    message( "logL improvement:", improvement )

  }   ## End of EM loop

  colnames(gamma) = tissues             ### assigning tissue names.
  results = list(gamma=gamma, alfa=alfa, beta=beta, sigma_g = sigma_g, sigma_e=sigma_e, m = m, logL=new_logL)

  message(Sys.time())
  message("-------MAPEM finished--------")

  return(results)

}


