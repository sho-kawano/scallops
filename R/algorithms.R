#' Optimize kernel parameters
#'
#' @param nburn Number of iterations
#' @param priors Prior distributions for gamma (Gaussian) and tau (Gamma)
#' @param s Coordinates of observations
#' @param y Values of observations
#' @param iso_kernel_matrix Isotropic kernel matrix
#'
#' @return List with optimized kernel parameters, as well as corresponding mean values of gamma and tau
#' @export 
#'
#' @examples
#' s = expand.grid(lat = c(-1:1), lon = c(-1:1))
#' y = c(1,0,0,0,1,0,0,0,1)
#' dpc_grid = get_grid(c(-1,1), c(-1,1), spacing = 2)
#' priors = get_priors(dpc_grid)
#' iso_kernel_matrix = get_kernel_matrix(s, dpc_grid)
#' fit = get_em(10, priors, s, y, iso_kernel_matrix)
#' cbind(get_estimates(s, dpc_grid, fit), obs = y)
#' 
get_em = function(nburn, priors, s, y, iso_kernel_matrix, dpc_grid) {
  s = get_mat2list(s)
  neighborsI = iso_kernel_matrix$neighborsI
  neighborsJ = iso_kernel_matrix$neighborsJ
  ngridpoints = iso_kernel_matrix$ngridpoints
  # initial values of rho, gamma, and tau
  rho = lapply(1:ngridpoints, FUN = function(i) c(0.25, 1, 1, 0))
  phi = get_phi(rho, iso_kernel_matrix)
  gamma = priors$gamma$mean
  tau = priors$tau$b / priors$tau$a
  kernel_matrix = get_kernel_matrix(s = s, dpc_grid = dpc_grid, phi = phi)
  # slight improvements on tau and gamma
  gamma = get_gamma_sample(priors$gamma$mean, priors$gamma$precision, 
                           kernel_matrix, tau, y, return_mean = TRUE)
  tau = get_tau_sample(kernel_matrix, gamma, y, priors$tau$a, priors$tau$b, return_mean = TRUE)
  # iterative method to optimize kernel parameters and obtain posterior means of gamma and tau
  cat('EM stage:\n')
  for (i in 1:nburn) {
    cat(i, ' ')
    for (j in 1:ngridpoints) {
      rho[[j]] = get_optimum_rho_j(j, dpc_grid, s, y, gamma, tau, rho, iso_kernel_matrix)
    }
    phi = get_phi(rho, iso_kernel_matrix)
    kernel_matrix = get_kernel_matrix(s = s, dpc_grid = dpc_grid, phi = phi)
    tau = get_tau_sample(kernel_matrix, gamma, y, priors$tau$a, priors$tau$b, return_mean = TRUE)
    gamma = get_gamma_sample(priors$gamma$mean, priors$gamma$precision, 
                             kernel_matrix, tau, y, return_mean = TRUE)
  }
  cat('\n')
  list(gamma = gamma, tau = tau, rho = rho)
}

#' Markov chain Monte Carlo (Gibbs) method to sample from the conditional distributions of gamma and tau
#'
#' @param s Observation locations
#' @param dpc_grid Discrete Process Convolution grid
#' @param y Observation values
#' @param nburn Number of EM iterations (burn-in) use dto estimate kernel parameters (rho)
#' @param nsample Number of MCMC iterations
#' @param priors Priors for gamma and tau
#' @param printEvery Console notifications every X iterations 
#' @param seed Seed for reproducible chains
#'
#' @return List with sampled values of gamma and tau, as well as optimal estimates for rho
#' @export
#'
#' @examples
#' s = data.frame(lat = rep(0, 3), lon = c(-1:1))
#' y = c(1,1,0)
#' dpc_grid = get_grid(c(-1,1), c(0,0), spacing = 2)
#' priors = get_priors(dpc_grid)
#' iso_kernel_matrix = get_kernel_matrix(s, dpc_grid)
#' fit = get_mcmc(s, dpc_grid, y, 10, 1000, priors, 100, 1)
#' cbind(get_estimates(s, dpc_grid, fit), obs = y)
get_mcmc = function(s, dpc_grid, y, nburn = 10, nsample = 10000, priors = NULL, printEvery = 1000, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  s = get_mat2list(s)
  iso_kernel_matrix = get_kernel_matrix(s = s, dpc_grid = dpc_grid)
  if (is.null(priors)) 
    priors = get_priors(dpc_grid)
  em = if (nburn > 0)
    get_em(nburn, priors, s, y, iso_kernel_matrix, dpc_grid)
  else
    list(gamma = priors$gamma$mean, 
         tau = priors$tau$b / priors$tau$a, 
         rho = lapply(1:iso_kernel_matrix$ngridpoints, FUN = function(i) c(0.25, 1, 1, 0)))
  if (nsample <= 0)
    em
  else {
    cat('MCMC stage:\n')
    gamma = em$gamma
    tau = em$tau
    phi = get_phi(em$rho, iso_kernel_matrix)
    kernel_matrix = get_kernel_matrix(s, dpc_grid, phi)
    sample_gamma = vector('list', length = nsample)
    sample_tau = vector('list', length = nsample)
    for (i in 1:nsample) {
      if (i %% printEvery == 0) cat(i, ' ')
      tau = get_tau_sample(kernel_matrix, gamma, y, priors$tau$a, priors$tau$b)
      gamma = get_gamma_sample(priors$gamma$mean, priors$gamma$precision, kernel_matrix, tau, y)
      sample_gamma[[i]] = gamma
      sample_tau[[i]] = tau
    }
    list(gamma = sample_gamma, tau = sample_tau, rho = em$rho)
  }
}