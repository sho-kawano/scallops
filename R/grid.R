#' Produce a regular Discrete Process Convolution grid.
#'
#' @param longitude_bounds Minimum and maximum longitude
#' @param latitude_bounds Minimum and maximum latitude
#' @param spacing Grid resolution, i.e., distance between adjacent gridpoints
#'
#' @return List with gridded coodinates, as well as metadata
#' @export
#'
#' @examples
#' get_grid(c(-1,1), c(2,4), 1)
get_grid = function(longitude_bounds, latitude_bounds, spacing) {
  lon = seq(longitude_bounds[1], longitude_bounds[2], by = spacing)
  lat = seq(latitude_bounds[1], latitude_bounds[2], by = spacing)
  list(lon = lon, lat = lat, spacing = spacing, coord = expand.grid(lat = lat,lon = lon))
}