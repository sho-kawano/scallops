get_loglik_rho_j = function(rho_candidate, my_phi, incomplete_phi, my_s, my_y, my_coord, my_gamma, my_keval, tau, grid_spacing){
  if (length(my_y) == 0)
    return(0)
  else {
    new_phi = mapply(my_phi, incomplete_phi, my_keval, FUN = function(lvec, inc_phi, keval) {
      inc_phi + keval * rho_candidate
    }, SIMPLIFY = FALSE)
    estimates = mapply(my_s, my_coord, my_gamma, new_phi, FUN = function(s_i, s_j, gamma_vec, phi_i) {
      keval = get_kernel_evaluation(s_i, s_j, grid_spacing, phi_i = phi_i)
      sum(keval * gamma_vec) / sum(keval)
    })
    loglik = -0.5 * tau * sum((my_y - estimates)**2)
    return(loglik)
  }
}

get_phi = function(rho, iso_kernel_matrix) {
  rho_mat = matrix(nrow = 4, unlist(rho))
  phi_mat = Matrix::tcrossprod(iso_kernel_matrix$mat, rho_mat)
  phi = lapply(1:iso_kernel_matrix$nobs, FUN = function(i) phi_mat[i, ])
  phi
}

get_optimum_rho_j = function(j, dpc_grid, s, y, gamma, tau, rho, iso_kernel_matrix){
  neighborsI = iso_kernel_matrix$neighborsI
  neighborsJ = iso_kernel_matrix$neighborsJ
  s_j = dpc_grid$coord[j, ] # coordinates of gridpoint j
  obs_idx = neighborsJ[[j]] # indices of obs in the neighborhood of j
  #my_kernel_sum = kernel_sum[obs_idx] # sum of kernel evaluations for each observation in the neighborhood of j
  my_s = lapply(obs_idx, FUN = function(i) s[[i]]) # coordinates of obs in the neighborhood of j
  my_y = y[obs_idx] #values of obs in the neighborhood of j
  my_neighborsI = lapply(obs_idx, FUN = function(i) neighborsI[[i]]) #gridpoints in the neighborhood of each obs
  my_pos = mapply(my_neighborsI, FUN = function(ngI) which(ngI == j)) #position of j in my_neighborsI
  my_coord = lapply(my_neighborsI, FUN = function(ngI) dpc_grid$coord[ngI, ])
  my_gamma = lapply(my_neighborsI, FUN = function(ngI) gamma[ngI])
  my_rho = lapply(my_neighborsI, FUN = function(ngI) lapply(ngI, FUN = function(jj) rho[[jj]]))
  iso_kernel_vec = mapply(obs_idx, my_neighborsI, SIMPLIFY = FALSE, FUN = function(ii, jj) {
    iso_kernel_matrix$mat[ii, jj]
  }) #list of isotropic evals
  my_keval = mapply(my_pos, iso_kernel_vec, FUN = function(pos, kv) kv[pos])
  my_phi = mapply(iso_kernel_vec, my_rho, FUN = function(kv, ph) {
    rho_mat = matrix(nrow = 4, unlist(ph))
    as.numeric(rho_mat %*% kv) # 4*1 phi vector
  }, SIMPLIFY = FALSE)
  incomplete_phi = mapply(my_phi, my_pos, iso_kernel_vec, FUN = function(ph, pos, kv) {
    ph - kv[pos] * rho[[j]]
  }, SIMPLIFY = FALSE)
  my_opt = optim(par = runif(4), fn = get_loglik_rho_j,
                 my_phi = my_phi, incomplete_phi = incomplete_phi, my_s = my_s, my_y = my_y,
                 my_coord = my_coord, my_gamma = my_gamma, my_keval = my_keval,
                 tau = tau, grid_spacing = dpc_grid$spacing,
                 method = 'L-BFGS-B', lower = 0, upper = 1, control = list(fnscale = -1))
  opt_rho = my_opt$par
  opt_rho
}