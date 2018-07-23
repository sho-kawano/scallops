#' Haversine distance between two points on the Earth's surface
#'
#' @param s_i Point i, identified by latitude and longitude
#' @param s_j Point j
#'
#' @return Distance between points i and j, in km 
#' @export
#'
#' @examples
#' d1 = haversine(s_i = list(lat = 0, lon = 10), s_j = list(lat = 0, lon = 11))
#' d2 = haversine(s_i = list(lat = 10, lon = 10), s_j = list(lat = 11, lon = 10))
#' abs(2 * pi * 6371 / 360 - d1) < 0.001
#' abs(2 * pi * 6371 / 360 - d2) < 0.001
#' 
haversine = function(s_i, s_j){
  earth_radius = 6371
  2 * earth_radius * asin(sqrt(sin((s_i$lat - s_j$lat) * pi / 360) ** 2 +
                                 cos(s_i$lat * pi / 180) * cos(s_j$lat * pi / 180) *
                                 sin((s_i$lon - s_j$lon) * pi / 360) ** 2))
}

#' Lemos and Sanso' (2012) Inverse Sigma matrix
#'
#' @param phi_i Kernel parameters
#' @param grid_spacing Distance between adjacent grid points (a.k.a. resolution)
#'
#' @return A 2x2 matrix used to compute the distance between two points
#' @export
#'
#' @examples
get_Sigma_inverse = function(phi_i, grid_spacing) {
  semimin = grid_spacing * (1 + phi_i[2])
  semimaj = semimin + phi_i[3] * (2 * grid_spacing - semimin)
  ang = pi * phi_i[4]
  Psi = 0.5 * c(1/semimin**2 + 1/semimaj**2, 1/semimin**2 - 1/semimaj**2)
  SigmaInv = matrix(nrow = 2, c(Psi[1] + Psi[2] * cos(2 * ang),
                                Psi[2] * sin(2 * ang),
                                Psi[2] * sin(2 * ang),
                                Psi[1] - Psi[2] * cos(2 * ang)))
  SigmaInv
}

#' Get kernel evaluation(s) centered at s_i and evaluated at one or many gridpoints s_j
#'
#' @param s_i Point where the kernel is centered at
#' @param s_j Grid points(s) where the kernel is evaluated at
#' @param grid_spacing DPC grid resolution
#' @param phi_i Kernel parameters at location s_i
#'
#' @return Kernel evaluation(s) at the gridpoint(s) s_j
#' @export
#'
#' @examples
#' get_kernel_evaluation(s_i = list(lat = 0, lon = 0), s_j = list(lat = 0, lon = 0), 
#'                       grid_spacing = 1, phi_i = c(0, 0, 0, 0)) == 1
#' get_kernel_evaluation(s_i = list(lat = 0, lon = 0), s_j = list(lat = 0, lon = 1), 
#'                       grid_spacing = 1, phi_i = c(0, 0, 0, 0)) == 0
#' get_kernel_evaluation(s_i = list(lat = 0, lon = 0), s_j = list(lat = 0, lon = 1), 
#'                       grid_spacing = 2, phi_i = c(0, 0, 0, 0)) == 0.75
#' 
get_kernel_evaluation = function(s_i, s_j, grid_spacing, phi_i){
  SigmaInv = get_Sigma_inverse(phi_i = phi_i, grid_spacing = grid_spacing)
  delta_lon = s_i$lon - s_j$lon
  delta_lat = s_i$lat - s_j$lat
  distance = delta_lon**2 * SigmaInv[1,1] + 2 * delta_lon * delta_lat * SigmaInv[1,2] + delta_lat**2 * SigmaInv[2,2]
  base = (1 - distance) * (distance < 1)
  evaluation = base ** (1 + 4 * phi_i[1])
  evaluation
}

#' Obtain the vector of kernel evaluations for all the gridpoints in the neighborhood of observation s_i
#'
#' @param s_i Observation location
#' @param neighbors_i List of gridpoints in the neighborhood of s_i
#' @param dpc_grid Discrete Process Convolution grid
#' @param phi_i Kernel parameters
#'
#' @return List with normalized vector of kernel evaluations, as well as sum of unnormalized evaluations
#' @export
#'
#' @examples
#' s_i = list(lat = 0, lon = 0)
#' dpc_grid = get_grid(c(-1, 1), c(0, 0), 1)
#' neighbors_i = c(1,2,3)
#' phi_i = c(0.5, 1, 0, 0)
#' get_kernel_vector(s_i, neighbors_i, dpc_grid, phi_i)
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

#' Identify the gridpoints in the neighborhood of location s_i
#'
#' @param s_i Observation location
#' @param dpc_grid Discrete Process Convolution grid
#'
#' @return Vector of indices that identify the grispoint neighbors of s_i
#' @export
#'
#' @examples
#' s_i = list(lat = 0, lon = 0)
#' dpc_grid = get_grid(c(0, 5), c(0, 5), 1)
#' get_neighbors_i(s_i, dpc_grid)
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

#' Identify the observation locations that are in the neighborhood of each gridpoint
#'
#' @param neighborsI List of gridpoint neighbors to each observation 
#' @param ngridpoints Number of gridpoints
#'
#' @return List of observation locations that are neighbors to each gridpoint
#' @export
#'
#' @examples
#' neighborsI = list(list(1, 2), list(2, 3))
#' get_neighbors_of_gridpoints(neighborsI, 3)
get_neighbors_of_gridpoints = function(neighborsI, ngridpoints) {
  idI = unlist(neighborsI)
  idJ = unlist(lapply(1:length(neighborsI), FUN = function(i) rep(i, length(neighborsI[[i]]))))
  q = cbind(idJ, idI)
  out = lapply(1:ngridpoints, FUN = function(j) q[q[,2] == j, 1])
  out
}

#' Construct a sparse kernel matrix
#'
#' @param s Observation locations
#' @param dpc_grid Discrete Process Convolution grid
#' @param phi Optional kernel parameters
#'
#' @return List with sparse matrix (nrow = number of observations, ncol = number of gridpoints), plus
#'         metadata
#' @export
#'
#' @examples
#' s = data.frame(lon = 1:3, lat = 1:3)
#' dpc_grid = get_grid(c(1, 3), c(1, 3), 1)
#' get_kernel_matrix(s, dpc_grid)$mat
#' 
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