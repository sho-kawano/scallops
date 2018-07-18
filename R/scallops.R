
get_grid = function(longitude_bounds, latitude_bounds, spacing) {
  lon = seq(longitude_bounds[1], longitude_bounds[2], by = spacing)
  lat = seq(latitude_bounds[1], latitude_bounds[2], by = spacing)
  list(lon = lon, lat = lat, spacing = spacing, coord = expand.grid(lat = lat,lon = lon))
}

get_neighbors_i = function(s_i, dpc_grid) {
  nlon = length(dpc_grid$lon)
  nlat = length(dpc_grid$lat)
  lo_min = max(1, sum(dpc_grid$lon < s_i$lon) - 1)
  lo_max = min(nlon, lo_min + 3)
  la_min = max(1, sum(dpc_grid$lat < s_i$lat) - 1)
  la_max = min(nlat, la_min + 3)
  indices = as.numeric(mapply(lo_min:lo_max, FUN = function(lo) nlat * (lo - 1) + la_min:la_max))
  indices
}

haversine = function(s_i, s_j){
  earth_radius = 6371
  2 * earth_radius * asin(sqrt(sin((s_i$lat - s_j$lat) * pi / 360) ** 2 +
                                 cos(s_i$lat * pi / 180) * cos(s_j$lat * pi / 180) *
                                 sin((s_i$lon - s_j$lon) * pi / 360) ** 2))
}

get_isotropic_kernel_evaluation = function(s_i, s_j, lambda_i){
  distance = haversine(s_i, s_j)
  ratio = distance / lambda_i$range_km
  base = (1 - ratio ** 2) * (ratio < 1)
  evaluation = base ** lambda_i$power
  evaluation
}

get_kernel_vector = function(s_i, neighbors_i, dpc_grid, lambda_i) {
  s_j = dpc_grid$coord[neighbors_i,]
  evals = get_isotropic_kernel_evaluation(s_i, s_j, lambda_i)
  sum_ev = sum(evals)
  if (sum_ev <= 0) {
    cat("Evaluations are all zero!\n")
    cat("s_i = ", s_i$lon, s_i$lat, "\n")
    cat("neighbors_i = ", neighbors_i, "\n")
    stop("Fatal error")
  }
  evals / sum_ev
}

get_kernel_matrix = function(s, neighbors, dpc_grid, lambda) {
  ii = as.numeric(unlist(mapply(1:length(neighbors), FUN = function(i) rep(i, length(neighbors[[i]])))))
  jj = as.numeric(unlist(neighbors))
  kv = as.numeric(unlist(mapply(s, neighbors, lambda, FUN = function(s_i, neighbors_i, lambda_i) {
    get_kernel_vector(s_i, neighbors_i, dpc_grid, lambda_i)
  })))
  sparseMatrix(i = ii, j = jj, x = kv, dims = c(length(s), nrow(dpc_grid$coord)))
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
          -nearest_neighbor_factor * diag_value
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
  prior_precision + tau * Matrix::crossprod(kernel_matrix)
}

get_posterior_mean = function(prior_mean, prior_precision, posterior_precision, kernel_matrix, tau, y) {
  d = prior_precision %*% prior_mean + tau * Matrix::crossprod(kernel_matrix, y)
  Matrix::solve(posterior_precision, d)
}

get_gamma_sample = function(prior_mean, prior_precision, kernel_matrix, tau, y) {
  posterior_precision = get_posterior_precision(prior_precision, kernel_matrix, tau)
  posterior_mean = get_posterior_mean(prior_mean, prior_precision, posterior_precision, kernel_matrix, tau, y)
  n = nrow(posterior_precision)
  R = Matrix::chol(posterior_precision)
  as.numeric(solve(R, rnorm(n)) + posterior_mean)
}

get_sum_squared_residuals = function(kernel_matrix, gamma, y) {
  sum((y - kernel_matrix %*% gamma) ** 2)
}

get_tau_sample = function(kernel_matrix, gamma, y, prior_tau_a, prior_tau_b) {
  sse = get_sum_squared_residuals(kernel_matrix, gamma, y)
  posterior_tau_a = prior_tau_a + 0.5 * length(y)
  posterior_tau_b = prior_tau_b + 0.5 * sse
  rgamma(1, shape = posterior_tau_a, rate = posterior_tau_b)
}

mcmc = function(nburn, nsample, prior_tau_a, prior_tau_b, prior_mean,
                prior_precision, kernel_matrix, y) {
  gamma = prior_mean
  tau = prior_tau_b / prior_tau_a
  for (i in 1:nburn) {
    tau = get_tau_sample(kernel_matrix, gamma, y, prior_tau_a, prior_tau_b)
    gamma = get_gamma_sample(prior_mean, prior_precision, kernel_matrix, tau, y)
  }
  sample_gamma = vector('list', length = nsample)
  sample_tau = vector('list', length = nsample)
  for (i in 1:nsample) {
    tau = get_tau_sample(kernel_matrix, gamma, y, prior_tau_a, prior_tau_b)
    gamma = get_gamma_sample(prior_mean, prior_precision, kernel_matrix, tau, y)
    sample_gamma[[i]] = gamma
    sample_tau[[i]] = tau
  }
  list(tau = sample_tau, gamma = sample_gamma)
}

test3 = function(){
  set.seed(123456)
  kernel_matrix = rbind(diag(1,5), diag(1,5), diag(1,5), diag(1,5), diag(1,5), diag(1,5))
  s = lapply(1:30, FUN = function(i) list(lat = 30 + runif(1, 0, 10), lon = -20 + runif(1, 0, 10)))
  y = as.numeric(mapply(s, FUN = function(s_i) {
    cos(s_i$lon * pi / 5) + cos(s_i$lat * pi / 5) + rnorm(1, 0, 0.01)
  }))
  dpc_grid = get_grid(c(-20, -10), c(30, 40), 2)
  ngrid = nrow(dpc_grid$coord)
  neighbors = lapply(s, FUN = function(s_i) get_neighbors_i(s_i, dpc_grid))
  lambda = lapply(s, FUN = function(s_i) list(range_km = 400, power = 2))
  kernel_matrix = get_kernel_matrix(s, neighbors, dpc_grid, lambda)
  nburn = 100
  nsample = 30
  prior_tau_a = 2
  prior_tau_b = 0.01
  prior_mean = rep(0, ngrid)
  prior_precision = diag(0.01, ngrid)
  mcmc_res = mcmc(nburn, nsample, prior_tau_a, prior_tau_b, prior_mean, prior_precision, kernel_matrix, y)
  plot(1 / sqrt(unlist(mcmc_res$tau)), ty = 'l', ylab = expression(paste("1 / ", sqrt(tau))))
  gm = Matrix(nrow = length(mcmc_res$gamma[[1]]), unlist(mcmc_res$gamma))
  image(gm)
}

test2 = function(){
  set.seed(123456)
  kernel_matrix = rbind(diag(1,5), diag(1,5), diag(1,5), diag(1,5), diag(1,5), diag(1,5))
  y = as.numeric(mapply(1:6, FUN = function(i) -2 + 0:4 + rnorm(5, sd = 0.01)))
  nburn = 1000
  nsample = 50
  prior_tau_a = 2
  prior_tau_b = 0.01
  prior_mean = rep(0, 5)
  prior_precision = diag(0.01, 5)
  mcmc_res = mcmc(nburn, nsample, prior_tau_a, prior_tau_b, prior_mean, prior_precision, kernel_matrix, y)
  plot(1 / sqrt(unlist(mcmc_res$tau)), ty = 'l', ylab = expression(paste("1 / ", sqrt(tau))))
  gm = Matrix(nrow = 5, unlist(mcmc_res$gamma))
  image(gm)
}


test = function(){
  set.seed(123456)
  dpc_grid = get_grid(c(-20, -10), c(30, 40), 2)
  n = 20
  s = lapply(1:n, FUN = function(x) list(lon = runif(1, -20, -10), lat = runif(1, 30, 40)))
  neighbors = lapply(s, FUN = function(s_i) get_neighbors_i(s_i, dpc_grid))
  lambda = lapply(1:n, FUN = function(i) list(range_km = 400, power = 2))
  km = get_kernel_matrix(s, neighbors, dpc_grid, lambda)
}

