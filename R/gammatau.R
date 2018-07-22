get_priors = function(dpc_grid, precision_diagonal_value = 1e-6, precision_nearest_neighbor_factor = 0) {
  ngridpoints = nrow(dpc_grid$coord)
  precision = if (precision_nearest_neighbor_factor == 0)
    Diagonal(ngridpoints, x = precision_diagonal_value)
  else
    get_prior_precision_matrix(dpc_grid, precision_diagonal_value, precision_nearest_neighbor_factor)
  list(gamma = list(mean = rep(0, ngridpoints), precision = precision), tau = list(a = 2, b = 2))
}

get_prior_precision_matrix = function(dpc_grid, diag_value, nearest_neighbor_factor) {
  nlon = length(dpc_grid$lon)
  nlat = length(dpc_grid$lat)
  ii = as.numeric(unlist(mapply(1:nlon, FUN = function(lo) {
    lo_min = max(1, lo - 1)
    lo_max = min(nlon, lo + 1)
    nx = lo_max - lo_min + 1
    unlist(mapply(1:nlat, FUN = function(la) {
      la_min = max(1, la - 1)
      la_max = min(nlat, la + 1)
      ny = la_max - la_min + 1
      idx = (lo - 1) * nlat + la
      rep(idx, nx * ny)
    }))
  })))
  jj = as.numeric(unlist(mapply(1:nlon, FUN = function(lo) {
    lo_min = max(1, lo - 1)
    lo_max = min(nlon, lo + 1)
    nx = lo_max - lo_min + 1
    unlist(mapply(1:nlat, FUN = function(la) {
      la_min = max(1, la - 1)
      la_max = min(nlat, la + 1)
      ny = la_max - la_min + 1
      idx = (lo - 1) * nlat + la
      as.numeric(mapply(lo_min:lo_max, FUN = function(llo) (llo - 1) * nlat + la_min:la_max))
    }))
  })))
  xx = as.numeric(unlist(mapply(1:nlon, FUN = function(lo) {
    lo_min = max(1, lo - 1)
    lo_max = min(nlon, lo + 1)
    nx = lo_max - lo_min + 1
    unlist(mapply(1:nlat, FUN = function(la) {
      la_min = max(1, la - 1)
      la_max = min(nlat, la + 1)
      ny = la_max - la_min + 1
      idx = (lo - 1) * nlat + la
      as.numeric(mapply(lo_min:lo_max, FUN = function(llo) mapply(la_min:la_max, FUN = function(lla) {
        if (llo == lo & lla == la)
          diag_value
        else if (llo != lo & lla != la)
          0
        else
          -abs(nearest_neighbor_factor) * diag_value
      })))
    }))
  })))
  mat = sparseMatrix(i = ii, j = jj, x = xx)
  mdet = Matrix::determinant(mat)
  if (mdet$sign == -1) {
    cat("Bad determinant, pick a smaller value for nearest_neighbor_factor.\n")
  } else {
    return(mat)
  }
}

get_posterior_precision = function(prior_precision, kernel_matrix, tau) {
  prior_precision + tau * Matrix::crossprod(kernel_matrix$mat)
}

get_posterior_mean = function(prior_mean, prior_precision, posterior_precision, kernel_matrix, tau, y) {
  d = prior_precision %*% prior_mean + tau * Matrix::crossprod(kernel_matrix$mat, y)
  Matrix::solve(posterior_precision, d)
}

get_gamma_sample = function(prior_mean, prior_precision, kernel_matrix, tau, y, return_mean = FALSE) {
  posterior_precision = get_posterior_precision(prior_precision, kernel_matrix, tau)
  posterior_mean = get_posterior_mean(prior_mean, prior_precision, posterior_precision, kernel_matrix, tau, y)
  n = nrow(posterior_precision)
  R = Matrix::chol(posterior_precision)
  if (return_mean)
    posterior_mean
  else
    as.numeric(solve(R, rnorm(n)) + posterior_mean)
}

get_tau_sample = function(kernel_matrix, gamma, y, prior_tau_a, prior_tau_b, return_mean = FALSE) {
  sse = sum((y - kernel_matrix$mat %*% gamma) ** 2)
  posterior_tau_a = prior_tau_a + 0.5 * length(y)
  posterior_tau_b = prior_tau_b + 0.5 * sse
  if (return_mean)
    posterior_tau_a / posterior_tau_b
  else
    rgamma(1, shape = posterior_tau_a, rate = posterior_tau_b)
}
