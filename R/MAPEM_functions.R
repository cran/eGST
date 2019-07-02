

## M-step: update the gamma matirx: for nTissues number of mixture model

optim_gamma = function(N, nTissues, Y, X, alfa, beta, w, sigma_e){

  gamma = matrix(0, N, nTissues)

  mean = as.list(1:nTissues)

  for(j in 1:nTissues) mean[[j]] = alfa[j] + as.numeric(X[[j]] %*% beta[[j]])

  log_comp = as.list(1:nTissues); total = numeric(N);
  log_w = log(w)

  for(j in 1:nTissues){
    log_comp[[j]] = log_w[j] + dnorm(Y, mean = mean[[j]], sd = rep(sigma_e[j], N), log = TRUE)
    total = total + exp(log_comp[[j]])
  }

  log_total = log(total); logL = (1/N)*sum(log_total);

  for(j in 1:nTissues) gamma[ ,j] = exp(log_comp[[j]] - log_total)

  results = list(gamma = gamma, logL = logL)

  return(results)

}


## M-step: update n (sum of gamma values across all indivs corresponding to a mixture component)

optim_n = function(gamma){
  n = colSums(gamma)
  return(n)
}

## M-step: weight for the sub-classes

optim_w = function(n){

  total_n = sum(n)
  w = n/total_n
  return(w)

}


## M-step: alfa vector across tissues.

optim_alfa = function(nTissues, N, n, gamma, Y, X, beta, sigma_e, sigma_alfa0){

  intercept = numeric(nTissues)
  total = 0

  for(j in 1:nTissues){

    total = sum( gamma[ ,j] * (Y - as.numeric(X[[j]] %*% beta[[j]])) )
    correct_factor = (sigma_e[j] / sigma_alfa0)^2
    intercept[j] = total / (n[j] + correct_factor)

  }

  return(intercept)

}


## M-step: beta (corresponding to j-th mixture component (tissue))

optim_beta = function(j, M, N, Geno, Gamma, Alfa, pheno, sigma_e, sigma_g){

  alfa = Alfa[j]                  ## Alfa: Vector of alfa params across tissues.
  gamma = Gamma[ ,j]              ## Gamma: matrix of posterior probabilities across tissues
  m = M[j]                        ## M: vector containing number of eQTLs specific to various tissues.
  X = Geno[[j]]                   ## Geno: List containing all genotype matrices across tissues.

  Y = pheno - alfa                ## pheno: phenotype vector

  ## creating the gamma0 matrix to multiply (element-wise multiplication) with X0 matrix and get the Z0 matrix

  sqrt_gamma = sqrt(gamma)
  sqrt_gamma = matrix(rep(sqrt_gamma, m), N, m)
  Z = sqrt_gamma*X
  rm(sqrt_gamma)
  total_matr = ( t(Z) %*% Z ) + ( ( (sigma_e[j]/sigma_g[j])^2 ) * diag(m) )
  rm(Z)

  W = matrix(rep((gamma*Y), m), N, m)
  V = W*X
  rm(W)
  total_vec = colSums(V)
  rm(V)

  #epsilon = 10^(-5); count=0;
  #while( (class(try(solve(total_matr),silent=T))=="matrix") == FALSE ){
  #  diag(total_matr) = diag(total_matr) + epsilon
  #  count = count+1;
  #}

  #inv = solve(total_matr)
  #Beta = inv %*% matrix(total_vec, m, 1)

  Beta = solve(total_matr, total_vec)

  return(Beta)
}


## M-step: sigma (error) vector across tissues

optim_sigma_e = function(nTissues, N, pheno, gamma, X, alfa, beta, n, a_e, b_e){

  Sigma = numeric(nTissues)

  for(j in 1:nTissues){

    Y = pheno - alfa[j]      ## pheno: phenotype vector
    total = 0

    diff = Y - as.numeric(X[[j]] %*% beta[[j]])
    total = sum(gamma[ ,j] * (diff^2))

    hyper_factor_numerator = 2*b_e
    hyper_factor_denominator = (2*a_e)+1

    numerator = total + hyper_factor_numerator
    denominator = n[j] + hyper_factor_denominator

    Sigma[j] = sqrt(numerator / denominator)

  }

  return(Sigma)

}


###-------- function to optimize for tau M-Step ----------###

tau_function = function(tau, m, si, beta, sigma_e){

  ((si-1)*log(tau)) - ( (m/2) * ( log(tau) + ( (as.numeric(t(beta)%*%beta)/(sigma_e^2)) * (1/tau) ) ) )

}


#### ============== M step: update tau parameters =============== ####

# optim_tau = function( nTissues, m, beta, sigma_e, si, epsilon ){
#
#   Tau = numeric(nTissues)
#
#   for(j in 1:nTissues){
#
#     res = optimize(f = tau_function, interval = c(epsilon, 1), maximum = TRUE, m = m[j], si = si, beta = beta[[j]], sigma_e = sigma_e[j] )
#
#     Tau[j] = res$maximum
#
#   }
#
#   return(Tau)
#
# }


###--------------- M step: update sigma_g ----------------###

optim_sigma_g = function( nTissues, m, beta, sigma_e, a_g, b_g ){

  Sigma_g = numeric(nTissues)

  for(j in 1:nTissues){

    nume = as.numeric( t(beta[[j]]) %*% beta[[j]] ) + (2*b_g)
    deno = m[j] + (2*a_g) + 1
    Sigma_g[j] = sqrt(nume/deno)

  }

  return(Sigma_g)

}



##---------- Compute the prior log-likelihood -------------##

log_prior_likelihood = function(nTissues, m, alfa, beta, sigma_alfa0, sigma_g, sigma_e, a_g, b_g, a_e, b_e){

  logf_alfa = numeric(nTissues)
  logf_beta = numeric(nTissues)
  logf_sigma_g = numeric(nTissues)
  logf_sigma_e = numeric(nTissues)

  for(j in 1:nTissues){

    logf_alfa[j] = - log(sigma_alfa0) - ( (1/(2*sigma_alfa0^2)) * alfa[j]^2 )
    logf_beta[j] = - (m[j]*log(sigma_g[j])) - ( (1/(2*sigma_g[j]^2)) * as.numeric( t(beta[[j]]) %*% beta[[j]] ) )
    logf_sigma_g[j] = - ( (2*a_g + 1)*log(sigma_g[j]) - (b_g/sigma_g[j]^2) )
    logf_sigma_e[j] = - ( (2*a_e + 1)*log(sigma_e[j]) - (b_e/sigma_e[j]^2) )

  }

  prior_logL = sum( logf_alfa + logf_beta + logf_sigma_g + logf_sigma_e )

  return( prior_logL )

}



#### compute the whole likelihood value for stopping EM iterations

full_log_likelihood = function(nTissues, w, N, X, Y, alfa, beta, sigma){

  mean = as.list(1:nTissues)

  for(j in 1:nTissues) mean[[j]] = alfa[j] + as.numeric(X[[j]] %*% beta[[j]])

  log_comp = as.list(1:nTissues); total = numeric(N);
  log_w = log(w)

  for(j in 1:nTissues){
    log_comp[[j]] = log_w[j] + dnorm(Y, mean = mean[[j]], sd = rep(sigma[j], N), log = TRUE)
    total = total + exp(log_comp[[j]])
  }

  log_total = log(total)

  logL = (1/N)*sum(log_total)

}





