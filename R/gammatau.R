#' Build basic priors for gamma and tau
#'
#' @param dpc_grid Discrete Process Convolution grid
#' @param precision_diagonal_value Diagonal value for gamma's prior precision matrix
#' @param precision_nearest_neighbor_factor Off-diagonal factor for gamma's prior precision matrix, in (0,1)
#'
#' @return List with prior mean and precision for gamma, and prior shape and rate parameters for tau
#' @export
#'
#' @examples
#' dpc_grid = get_grid(c(0,1), c(1,2), 1)
#' priors = get_priors(dpc_grid, 5, 0.4)
#' priors$gamma$precision
#' (prior_variance = solve(priors$gamma$precision))
#' 
get_priors = function(dpc_grid, precision_diagonal_value = 1e-6, 
                      precision_nearest_neighbor_factor = 0) {
  ngridpoints = nrow(dpc_grid$coord)
  precision = if (precision_nearest_neighbor_factor == 0)
    Diagonal(ngridpoints, x = precision_diagonal_value)
  else
    get_prior_precision_matrix(dpc_grid, precision_diagonal_value, precision_nearest_neighbor_factor)
  list(gamma = list(mean = rep(0, ngridpoints), precision = precision), tau = list(a = 2, b = 2))
}

#' Construct gamma's precision matrix for a given DPC grid
#'
#' @param dpc_grid Discrete Process Convolution grid
#' @param diag_value Diagonal value in precision matrix
#' @param nearest_neighbor_factor Off-diagonal factor for adjacent gridpoints
#'
#' @return A sparse precision matrix
#'
#' @examples
#' dpc_grid = get_grid(c(0,1), c(1,2), 1)
#' get_prior_precision_matrix(dpc_grid, 3, 0.2)
#' 
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

#' Compute the full conditional (given tau, rho, and data) precision of gamma
#'
#' @param prior_precision Prior precision matrix for gamma
#' @param kernel_matrix Discrete Process Convolution kernel matrix (fixed values of rho)
#' @param tau Sampled value of tau
#'
#' @return Full conditional precision of gamma
#' @export
#'
#' @examples
#' get_posterior_precision(Matrix::Diagonal(2, 3), list(mat = Matrix::Diagonal(2, 1)), 5)
#' 
get_posterior_precision = function(prior_precision, kernel_matrix, tau) {
  prior_precision + tau * Matrix::crossprod(kernel_matrix$mat)
}

#' Compute the full conditional (given tau, rho, and data) mean of gamma
#'
#' @param prior_mean Prior mean of gamma
#' @param prior_precision Prior precision of gamma
#' @param posterior_precision Full conditional precision of gamma
#' @param kernel_matrix DPC kernel matrix
#' @param tau Sampled value of tau
#' @param y Data
#'
#' @return
#' @export
#'
#' @examples
#' km = list(mat = Matrix::Diagonal(2, 1))
#' prior_pr = Matrix::Diagonal(2, 0.3)
#' tau = 5
#' post_pr = get_posterior_precision(prior_pr, km, tau)
#' get_posterior_mean(prior_mean = rep(0, 2), prior_precision = prior_pr, 
#'                    posterior_precision = post_pr, kernel_matrix = km, tau = tau, y = c(0.9, 0.7))
get_posterior_mean = function(prior_mean, prior_precision, posterior_precision, kernel_matrix, tau, y) {
  d = prior_precision %*% prior_mean + tau * Matrix::crossprod(kernel_matrix$mat, y)
  Matrix::solve(posterior_precision, d)
}

#' Sample gamma from its full conditional distribution
#'
#' @param prior_mean Prior mean of gamma
#' @param prior_precision Prior precision of gamma
#' @param kernel_matrix DPC kernel matrix
#' @param tau Sampled value of tau
#' @param y Observation value(s)
#' @param return_mean Provide full conditional mean as pseudo-sample?
#'
#' @return Either the full conditional mean of gamma or a sample of it
#' @export
#'
#' @examples
#' km = list(mat = Matrix::Diagonal(2, 1))
#' prior_pr = Matrix::Diagonal(2, 0.3)
#' tau = 5
#' post_pr = get_posterior_precision(prior_pr, km, tau)
#' get_gamma_sample(prior_mean = rep(0, 2), prior_precision = prior_pr, 
#'                  kernel_matrix = km, tau = tau, y = c(0.9, 0.7), return_mean = TRUE)
#' 
get_gamma_sample = function(prior_mean, prior_precision, kernel_matrix, tau, y, return_mean = FALSE) {
  posterior_precision = get_posterior_precision(prior_precision, kernel_matrix, tau)
  posterior_mean = get_posterior_mean(prior_mean, prior_precision, posterior_precision, kernel_matrix, tau, y)
  n = nrow(posterior_precision)
  if (return_mean)
    as.numeric(posterior_mean)
  else {
    R = Matrix::chol(posterior_precision)
    as.numeric(solve(R, rnorm(n)) + posterior_mean)
  }
}

#' Sample tau from its full conditional distribution
#'
#' @param kernel_matrix DPC kernel matrix
#' @param gamma Sample of gamma
#' @param y Observation values
#' @param prior_tau_a Prior shape parameter
#' @param prior_tau_b Prior rate parameter
#' @param return_mean Provide full conditional mean as pseudo-sample?
#'
#' @return Either the full conditional mean of tau or a sample of it
#' @export
#'
#' @examples
#' km = list(mat = Matrix::Diagonal(2, 1))
#' get_tau_sample(kernel_matrix = km, gamma = c(0.849, 0.66), y = c(0.9, 0.7), 
#'                prior_tau_a = 2, prior_tau_b = 2, return_mean = TRUE)
#'                  
get_tau_sample = function(kernel_matrix, gamma, y, prior_tau_a, prior_tau_b, return_mean = FALSE) {
  sse = sum((y - kernel_matrix$mat %*% gamma) ** 2)
  posterior_tau_a = prior_tau_a + 0.5 * length(y)
  posterior_tau_b = prior_tau_b + 0.5 * sse
  if (return_mean)
    posterior_tau_a / posterior_tau_b
  else
    rgamma(1, shape = posterior_tau_a, rate = posterior_tau_b)
}
