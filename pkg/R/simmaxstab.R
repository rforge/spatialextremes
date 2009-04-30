rmaxstab <- function(n, coord, model = "gauss", cov11 = 1, cov12 = 0,
                     cov22 = 1, grid = FALSE, ...){

  if (model != "gauss")
    stop("Currently this function is only available for the Smith model")
  
  dist.dim <- ncol(coord)

  if (dist.dim != 2)
    stop("Currently this function is only available for R^2")
  
  n.site <- nrow(coord)
  coord.range <- apply(coord, 2, range)
  center <- colMeans(coord.range)
  edge <- max(apply(coord.range, 2, diff))
   
  if (grid){
    new.coord <- NULL
    n.grid <- nrow(coord)
    n.site <- n.grid^2
    n <- 1

    for (x in coord[,1])
      new.coord <- rbind(new.coord, cbind(x, coord[,2]))

    coord <- new.coord
    n.site <- nrow(coord)
  }

  ans <- .C("rsmith", as.double(coord), as.double(center), as.double(edge),
            as.integer(n), as.integer(n.site), as.double(cov11), as.double(cov12),
            as.double(cov22), ans = double(n.site * n), PACKAGE = "SpatialExtremes")$ans

  if (grid)
    ans <- matrix(ans, n.grid, byrow = TRUE)

  else
    ans <- matrix(ans, nrow = n, ncol = n.site, byrow = TRUE)
  
  return(ans)
}
