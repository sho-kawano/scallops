get_grid = function(longitude_bounds, latitude_bounds, spacing) {
  lon = seq(longitude_bounds[1], longitude_bounds[2], by = spacing)
  lat = seq(latitude_bounds[1], latitude_bounds[2], by = spacing)
  list(lon = lon, lat = lat, spacing = spacing, coord = expand.grid(lat = lat,lon = lon))
}