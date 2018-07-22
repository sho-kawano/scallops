get_em = function(nburn, priors, s, y, iso_kernel_matrix) {
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


get_mcmc = function(s, dpc_grid, y, nburn = 10, nsample = 10000, priors = NULL, printEvery = 1000, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  s = get_mat2list(s)
  iso_kernel_matrix = get_kernel_matrix(s = s, dpc_grid = dpc_grid)
  if (is.null(priors)) 
    priors = get_priors(dpc_grid)
  em = if (nburn > 0)
    get_em(nburn, priors, s, y, iso_kernel_matrix)
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
    list(tau = sample_tau, gamma = sample_gamma, rho = em$rho)
  }
}