get_gridpoint_influence = function(dpc_grid, lat, lon, fit = NULL) {
  idx = which(dpc_grid$coord$lat == lat & dpc_grid$coord$lon == lon)
  if (length(idx) < 1) {
    cat('Could not find gridpoint.\n')
    return()
  }
  s = expand.grid(lat = seq(dpc_grid$lat[1], dpc_grid$lat[length(dpc_grid$lat)], by = dpc_grid$spacing / 8),
                  lon = seq(dpc_grid$lon[1], dpc_grid$lon[length(dpc_grid$lon)], by = dpc_grid$spacing / 8))
  kernel_par = if (is.null(fit)) c(0.25, 1, 1, 0) else fit$rho[[idx]]
  phi = lapply(1:nrow(s), FUN = function(i) kernel_par)
  iso_kernel_matrix = get_kernel_matrix(s = s, dpc_grid = dpc_grid, phi = phi)
  z = as.numeric(iso_kernel_matrix$mat[, idx])
  df = cbind(s, z)
  df
}

#' Obtain dataframe of kernel ellipses centered at the DPC gridpoints
#'
#' @param dpc_grid Discrete Process Convolution grid
#' @param fit Fitted model
#' @param eval Size of ellipses: k[s_i, s_j, \phi] = eval
#'
#' @return Dataframe with coordinates of points that are points in the ellipses; gp identifies the gridpoint
#'
#' @examples
get_ellipses = function(dpc_grid, fit, eval) {
  angles = seq(0, 2 * pi, length = 100)
  ngridpoints = nrow(dpc_grid$coord)
  ell = mapply(1:ngridpoints, FUN = function(j) {
    center = dpc_grid$coord[j, ]
    rho = fit$rho[[j]]
    SigmaInv = get_Sigma_inverse(phi_i = rho, grid_spacing = dpc_grid$spacing)
    aux = cos(angles)**2 * SigmaInv[1,1] + 2 * sin(angles) * cos(angles) * SigmaInv[1,2] + 
      sin(angles)**2 * SigmaInv[2,2]
    pow = 1 + 4 * rho[1]
    r = (1 - eval^(1/pow)) / sqrt(aux)
    cbind(center$lon + r * cos(angles), center$lat + r * sin(angles))
  })
  lon = as.numeric(ell[1:100,])
  lat = as.numeric(ell[101:200,])
  gp = as.numeric(mapply(1:ngridpoints, FUN = function(i) rep(i, 100)))
  out = data.frame(lon = lon, lat = lat, gp = gp)
  out
}

#' Get posterior predictive mean estimates at location(s) in the domain
#'
#' @param s Locations
#' @param dpc_grid Discrete Process Convolution grid
#' @param fit Fitted model
#'
#' @return Data frame with coordinates and posterior predictive means
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
#' 
get_estimates = function(s, dpc_grid, fit) {
  slist = get_mat2list(s)
  gamma_mat = if (class(fit$gamma) == 'dgeMatrix')
    as.matrix(fit$gamma)
  else
    matrix(nrow = length(fit$gamma[[1]]), unlist(fit$gamma))
  gamma_mean = apply(gamma_mat, MARGIN = 1, FUN = mean)
  iso_kernel_matrix = get_kernel_matrix(s = slist, dpc_grid = dpc_grid)
  rho_mat = matrix(nrow = 4, unlist(fit$rho))
  phi_mat = Matrix::tcrossprod(iso_kernel_matrix$mat, rho_mat)
  phi = lapply(1:nrow(phi_mat), FUN = function(i) phi_mat[i, ])
  km = get_kernel_matrix(s = slist, dpc_grid = dpc_grid, phi = phi)
  estimates = as.numeric(km$mat %*% gamma_mean)
  df = cbind(s, z = estimates)
  df
}

#' Get a point estimate for a single location
#'
#' @param lon Longitude of point
#' @param lat Latitude of point
#' @param dpc_grid Discrete Process Convolution grid
#' @param fit Fitted model
#'
#' @return Posterior predictive mean at given location
#' @export
#'
#' @examples
#' s = data.frame(lat = rep(0, 3), lon = c(-1:1))
#' y = c(1,1,0)
#' dpc_grid = get_grid(c(-1,1), c(0,0), spacing = 2)
#' priors = get_priors(dpc_grid)
#' iso_kernel_matrix = get_kernel_matrix(s, dpc_grid)
#' fit = get_mcmc(s, dpc_grid, y, 10, 1000, priors, 100, 1)
#' get_point_estimate(lon = s$lon[1], lat = s$lat[1], dpc_grid = dpc_grid, fit = fit)
#' 
get_point_estimate = function(lon, lat, dpc_grid, fit) {
  s = data.frame(lon = lon, lat = lat)
  get_estimates(s = s, dpc_grid = dpc_grid, fit = fit)
}

#' Plot gridded posterior predictive means
#'
#' @param obs_coord Coordinates of original observation locations
#' @param dpc_grid Discrete Process Convolution grid
#' @param fit Fitted model
#' @param fineness Resolution of gridded means (1 = coarse, 10 = fine)
#'
#' @return Gridded estimates within the domain defined by the observation locations (extapolation is dangerous)
#' @export
#'
#' @examples
#' set.seed(5)
#' s = data.frame(lat = runif(3, min = 0, max = 1), lon = runif(3, -1, 1))
#' y = c(1,1,0)
#' dpc_grid = get_grid(c(-1,1), c(-1,1), spacing = 2)
#' priors = get_priors(dpc_grid)
#' iso_kernel_matrix = get_kernel_matrix(s, dpc_grid)
#' fit = get_mcmc(s, dpc_grid, y, 10, 1000, priors, 100, 1)
#' gr = get_gridded_estimates(obs_coord = s, dpc_grid = dpc_grid, fit = fit, fineness = 10)
#' plot(gr$lon, gr$lat, cex = 0.1 + 4 * gr$z, xlim = c(-1, 1), ylim = c(0, 1))
#' points(s$lon, s$lat, pch = 19, cex = 0.1 + 4 * y, col = 'red')
#' 
get_gridded_estimates = function(obs_coord, dpc_grid, fit, fineness) {
  obs_mat = get_list2mat(obs_coord)
  triangulation = tri.mesh(obs_mat$lon, obs_mat$lat, duplicate="remove")
  s = expand.grid(lat = seq(dpc_grid$lat[1], dpc_grid$lat[length(dpc_grid$lat)], by = dpc_grid$spacing / fineness),
                  lon = seq(dpc_grid$lon[1], dpc_grid$lon[length(dpc_grid$lon)], by = dpc_grid$spacing / fineness))
  valid = in.convex.hull(tri.obj = triangulation, x = s$lon, y = s$lat)
  s = s[valid, ]
  get_estimates(s = s, dpc_grid = dpc_grid, fit = fit)
}


#' Plot the influence of a gridded Gaussian variable on the interpolation
#'
#' @param dpc_grid Discrete Process Convolution grid
#' @param lat Latitute of desired gridpoint
#' @param lon Longitude of desired gridpoint
#' @param fit (Optional) fitted model
#'
#' @return A ggplot showing the weight of the selected variable on the fit
#' @export
#'
#' @examples
#' get_influence_plot(get_grid(c(-1,1), c(-1,1), spacing = 2), -1, -1, fit = NULL)
get_influence_plot = function(dpc_grid, lat, lon, fit = NULL) {
  df = get_gridpoint_influence(dpc_grid, lat, lon, fit)
  p = ggplot() + geom_raster(data = df[df$z > 0,], aes(x = lon, y = lat, fill = z))  +
    geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red') +
    coord_fixed(ratio = 1) +
    scale_fill_gradient(low = "light gray", high = "dark blue", space = "Lab",
                        na.value = "grey50", guide = "colourbar", aesthetics = "fill")
  return(p)
}

#' Plot kernel ellipses centered at the DPC gridpoints
#'
#' @param dpc_grid Discrete Process Convolution grid
#' @param fit Fitted model
#' @param eval Size of ellipses: k[s_i, s_j, \phi] = eval
#'
#' @return Plot of ellipses
#' @export
#'
#' @examples
#' s = expand.grid(lat = -1:1, lon=-1:1)
#' y = c(1,0,0,0,1,0,0,0,1)
#' dpc_grid = get_grid(c(-1,1), c(-1,1), spacing = 2)
#' priors = get_priors(dpc_grid)
#' iso_kernel_matrix = get_kernel_matrix(s, dpc_grid)
#' fit = get_mcmc(s, dpc_grid, y, 10, 1000, priors, 100, 1)
#' get_interpolation_plot(s, dpc_grid, fit)
#' get_ellipses_plot(dpc_grid, fit)
get_ellipses_plot = function(dpc_grid, fit, eval = 0.7) {
  el = get_ellipses(dpc_grid, fit, eval)
  p = ggplot() + geom_path(data = el, aes(x = lon, y = lat, group = gp), color = 'red') +
    geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red') +
    coord_fixed(ratio = 1)
  return(p)
}

#' Plot the interpolated surface
#'
#' @param obs_coord Coordinates of original observations
#' @param dpc_grid Discrete Process Convolution grid
#' @param fit Fitted model
#' @param fineness Resolution of interpolation (1 = coarse, 10 = fine)
#' @param contour_binwidth Width of bins for contour plot
#' @param do_grid Plot DPC gridpoints?
#'
#' @return Plot of gridded posterior predictive means
#' @export
#'
#' @examples
#' s = expand.grid(lat = -1:1, lon=-1:1)
#' y = c(1,0,0,0,1,0,0,0,1)
#' dpc_grid = get_grid(c(-1,1), c(-1,1), spacing = 2)
#' priors = get_priors(dpc_grid)
#' iso_kernel_matrix = get_kernel_matrix(s, dpc_grid)
#' fit = get_mcmc(s, dpc_grid, y, 10, 1000, priors, 100, 1)
#' get_interpolation_plot(s, dpc_grid, fit, contour_binwidth = 0.1)
get_interpolation_plot = function(obs_coord, dpc_grid, fit, fineness = 16,
                                  contour_binwidth = NULL, do_grid = TRUE) {
  df = get_gridded_estimates(obs_coord = obs_coord, dpc_grid = dpc_grid, fit = fit, fineness)
  ctr = if (is.null(contour_binwidth))
    NULL
  else {
    geom_contour(data = df, aes(x = lon, y = lat, z = z), binwidth = contour_binwidth)
  }
  grd = if (do_grid)
    geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red')
  else
    NULL
  p = ggplot(data = df, mapping = aes(x = lon, y = lat)) + 
    geom_raster(data = df, mapping = aes(x = lon, y = lat, fill = z)) + 
    scale_fill_gradient2(low = "black", mid = 'light blue', high = "white", space = "Lab", 
                         midpoint = mean(df$z), na.value = "black", 
                         guide = "colourbar", aesthetics = "fill") +
    ctr + grd + coord_fixed(ratio = 1)
  return(p)
}