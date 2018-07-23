get_gridpoint_influence = function(dpc_grid, lat, lon, fit = NULL) {
  idx = which(dpc_grid$coord$lat == lat & dpc_grid$coord$lon == lon)
  if (length(idx) != 1) {
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

get_ellipses = function(dpc_grid, fit, eval) {
  angles = seq(0, 2 * pi, length = 100)
  ngridpoints = nrow(dpc_grid$coord)
  ell = mapply(1:ngridpoints, FUN = function(j) {
    center = dpc_grid$coord[j, ]
    rho = fit$rho[[j]]
    Sigma = get_Sigma(phi_i = rho, grid_spacing = dpc_grid$spacing)
    aux = cos(angles)**2 * Sigma[1,1] + 2 * sin(angles) * cos(angles) * Sigma[1,2] + sin(angles)**2 * Sigma[2,2]
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

get_point_estimate = function(lon, lat, dpc_grid, fit) {
  s = data.frame(lon = lon, lat = lat)
  get_estimates(s = s, dpc_grid = dpc_grid, fit = fit)
}

get_gridded_estimates = function(obs_coord, dpc_grid, fit, fineness) {
  obs_mat = get_list2mat(obs_coord)
  triangulation = tri.mesh(obs_mat$lon, obs_mat$lat, duplicate="remove")
  s = expand.grid(lat = seq(dpc_grid$lat[1], dpc_grid$lat[length(dpc_grid$lat)], by = dpc_grid$spacing / fineness),
                  lon = seq(dpc_grid$lon[1], dpc_grid$lon[length(dpc_grid$lon)], by = dpc_grid$spacing / fineness))
  valid = in.convex.hull(tri.obj = triangulation, x = s$lon, y = s$lat)
  s = s[valid, ]
  get_estimates(s = s, dpc_grid = dpc_grid, fit = fit)
}

get_influence_plot = function(dpc_grid, lat, lon, fit = NULL) {
  df = get_gridpoint_influence(dpc_grid, lat, lon, fit)
  ggplot() + geom_raster(data = df[df$z > 0,], aes(x = lon, y = lat, fill = z))  +
    geom_point(data = scallop, aes(x = longitude, y = latitude, size = logcatch), shape = 1) +
    geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red') +
    coord_fixed(ratio = 1) +
    scale_fill_gradient(low = "light gray", high = "dark blue", space = "Lab",
                        na.value = "grey50", guide = "colourbar", aesthetics = "fill")
}

get_ellipses_plot = function(dpc_grid, fit, eval = 0.7) {
  el = get_ellipses(dpc_grid, fit, eval)
  ggplot() + geom_path(data = el, aes(x = lon, y = lat, group = gp), color = 'red') +
    geom_point(data = scallop, aes(x = longitude, y = latitude, size = logcatch), shape = 1) +
    geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red') +
    coord_fixed(ratio = 1)
}

get_interpolation_plot = function(obs_coord, dpc_grid, fit, fineness = 16, obs_as_labels = FALSE,
                                  contour_binwidth = NULL, do_grid = TRUE) {
  df = get_gridded_estimates(obs_coord = scallop, dpc_grid = dpc_grid, fit = fit, fineness)
  pts = if (obs_as_labels)
    geom_text(data = scallop, aes(x = longitude, y = latitude, label = round(logcatch)))
  else
    geom_point(data = scallop, aes(x = longitude, y = latitude, size = logcatch), shape = 1)
  ctr = if (is.null(contour_binwidth))
    NULL
  else {
    geom_contour(data = df, aes(x = lon, y = lat, z = z), binwidth = contour_binwidth)
  }
  grd = if (do_grid)
    geom_point(data = dpc_grid$coord, aes(x = lon, y = lat), shape = 3, color = 'red')
  else
    NULL
  ggplot(data = df, mapping = aes(x = lon, y = lat)) + 
    geom_raster(data = df, mapping = aes(x = lon, y = lat, fill = z)) + 
    scale_fill_gradient2(low = "black", mid = 'light blue', high = "white", space = "Lab", 
                         midpoint = mean(df$z), na.value = "black", 
                         guide = "colourbar", aesthetics = "fill") +
    ctr + pts + grd + coord_fixed(ratio = 1)
}