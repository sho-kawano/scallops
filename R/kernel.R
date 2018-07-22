haversine = function(s_i, s_j){
  earth_radius = 6371
  2 * earth_radius * asin(sqrt(sin((s_i$lat - s_j$lat) * pi / 360) ** 2 +
                                 cos(s_i$lat * pi / 180) * cos(s_j$lat * pi / 180) *
                                 sin((s_i$lon - s_j$lon) * pi / 360) ** 2))
}

get_Sigma = function(phi_i, grid_spacing) {
  semimin = grid_spacing * (1 + phi_i[2])
  semimaj = semimin + phi_i[3] * (2 * grid_spacing - semimin)
  ang = pi * phi_i[4]
  Psi = 0.5 * c(1/semimin**2 + 1/semimaj**2, 1/semimin**2 - 1/semimaj**2)
  Sigma = matrix(nrow = 2, c(Psi[1] + Psi[2] * cos(2 * ang),
                             Psi[2] * sin(2 * ang),
                             Psi[2] * sin(2 * ang),
                             Psi[1] - Psi[2] * cos(2 * ang)))
  Sigma
}

get_kernel_evaluation = function(s_i, s_j, grid_spacing, phi_i){
  Sigma = get_Sigma(phi_i = phi_i, grid_spacing = grid_spacing)
  delta_lon = s_i$lon - s_j$lon
  delta_lat = s_i$lat - s_j$lat
  distance = delta_lon**2 * Sigma[1,1] + 2 * delta_lon * delta_lat * Sigma[1,2] + delta_lat**2 * Sigma[2,2]
  base = (1 - distance) * (distance < 1)
  evaluation = base ** (1 + 4 * phi_i[1])
  evaluation
}

get_kernel_vector = function(s_i, neighbors_i, dpc_grid, phi_i) {
  s_j = dpc_grid$coord[neighbors_i,]
  evals = get_kernel_evaluation(s_i, s_j, dpc_grid$spacing, phi_i)
  sum_ev = sum(evals)
  if (sum_ev <= 0) {
    cat("Evaluations are all zero!\n")
    cat("s_i = ", s_i$lon, s_i$lat, "\n")
    cat("neighbors_i = ", neighbors_i, "\n")
    stop("Fatal error")
  }
  return(list(normalized_evals = evals / sum_ev, eval_sum = sum_ev))
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

get_neighbors_of_gridpoints = function(neighborsI, ngridpoints) {
  idI = unlist(neighborsI)
  idJ = unlist(lapply(1:length(neighborsI), FUN = function(i) rep(i, length(neighborsI[[i]]))))
  q = cbind(idJ, idI)
  out = lapply(1:ngridpoints, FUN = function(j) q[q[,2] == j, 1])
  out
}

get_kernel_matrix = function(s, dpc_grid, phi = NULL) {
  s = get_mat2list(s)
  neighborsI = lapply(s, FUN = function(s_i) get_neighbors_i(s_i, dpc_grid))
  neighborsJ = get_neighbors_of_gridpoints(neighborsI, nrow(dpc_grid$coord))
  if (is.null(phi)) phi = lapply(1:length(s), FUN = function(i) c(0.25, 1, 1, 0))
  ii = as.numeric(unlist(mapply(1:length(neighborsI), FUN = function(i) rep(i, length(neighborsI[[i]])))))
  jj = as.numeric(unlist(neighborsI))
  kv = mapply(s, neighborsI, phi, FUN = function(s_i, neighbors_i, phi_i) {
    get_kernel_vector(s_i, neighbors_i, dpc_grid, phi_i)
  }, SIMPLIFY = FALSE)
  xx = as.numeric(unlist(lapply(kv, FUN = function(x) x$normalized_evals)))
  mat = sparseMatrix(i = ii, j = jj, x = xx, dims = c(length(s), nrow(dpc_grid$coord)))
  list(mat = mat, ngridpoints = ncol(mat), nobs = nrow(mat), neighborsI = neighborsI, neighborsJ = neighborsJ)
}