distance <- function(coord, vec = FALSE){
  ##This function computes the distance between each pair of locations

  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- nrow(coord)
  dist.dim <- ncol(coord)
  n.pairs <- n.site * (n.site - 1) / 2

  if (vec){
    dist <- .C("distVecFct", as.double(coord), as.integer(n.site),
               as.integer(dist.dim), distVec = double(dist.dim * n.pairs),
               PACKAGE = "SpatialExtremes")$distVec
    dist <- matrix(dist, ncol = dist.dim, nrow = n.pairs)
  }

  else    
    dist <- .C("distance", as.double(coord), as.integer(dist.dim),
               as.integer(n.site), dist = double(n.pairs),
               PACKAGE = "SpatialExtremes")$dist
  
  return(dist)
}

gev2frech <- function(x, loc, scale, shape)
  pmax(1 + shape * (x - loc) / scale, 0)^(1/shape)
