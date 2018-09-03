#' Convert matrix/data.frame observation locations into list
#'
#' @param s Observation locations
#'
#' @return List of obs locations
#'
get_mat2list = function(s) {
  if (is.data.frame(s) | is.matrix(s))
    lapply(1:nrow(s), FUN = function(i) s[i,])
  else
    s
}

#' Convert list of observation locations into data frame
#'
#' @param s List of observation locations
#'
#' @return Data frame of obs locations
#'
get_list2mat = function(s) {
  if (is.data.frame(s) | is.matrix(s))
    s
  else
    data.frame(lon = mapply(s, FUN = function(ss) ss$lon),
               lat = mapply(s, FUN = function(ss) ss$lat))
}