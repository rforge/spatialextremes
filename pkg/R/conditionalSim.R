condrgp <- function(n, coord, data.coord, data, cov.mod = "powexp",
                    mean = 0, nugget = 0, sill = 1, range = 1,
                    smooth = 1, grid = FALSE, control = list()){

  if (cov.mod == "caugen")
    stop("The generalized Cauchy covariance family isn't implemented")
  
  if (is.null(coord) & !is.null(data.coord))
    stop("'coord' and 'data.coord' don't match")

  if (!is.null(dim(coord)) & (any(ncol(coord) != ncol(data.coord))))
    stop("'coord' and 'data.coord' don't match")

  new.coord <- coord

  if (grid){
    if (is.null(dim(coord)))
      stop("You cannot use 'grid = TRUE' for 1 dimensional processes")
    
    dummy <- NULL
    for (i in 1:nrow(new.coord))
      dummy <- rbind(dummy, cbind(new.coord[,1], new.coord[i,2]))

    new.coord <- dummy
  }
  
  if (is.null(dim(coord))){
    n.cond <- length(data.coord)
    n.loc <- length(new.coord)
    new.coord <- c(data.coord, new.coord)
  }

  else{
    n.cond <- nrow(data.coord)
    n.loc <- nrow(coord)
    new.coord <- rbind(data.coord, new.coord)
  }
    
  uncond <- rgp(n, new.coord, cov.mod = cov.mod, mean = mean, nugget = nugget,
                sill = sill, range = range, smooth = smooth, control = control)

  weights <- kriging(data, data.coord, coord, cov.mod = cov.mod, sill = sill,
                     range = range, smooth = smooth, grid = grid,
                     only.weights = TRUE)$weights
  if (grid) {
    ans <- array(NA, c(n.loc, n.loc, n))

    for (i in 1:n){
      res <- data - uncond[i, 1:n.cond]
      krig <- matrix(res %*% weights, n.loc)
      ans[,,i] <- uncond[i,-(1:n.cond)] + krig
    }
  }

  else {
    ans <- matrix(NA, n, n.loc)
    
    for (i in 1:n){
      res <- data - uncond[i,1:n.cond]
      krig <- res %*% weights      
      ans[i,] <- uncond[i,-(1:n.cond)] + krig
    }
  }

  if (grid & (n == 1))
    ans <- matrix(ans, n.loc)

  return(list(coord = coord, cond.sim = ans, data.coord = data.coord,
              data = data, cov.mod = cov.mod, grid = grid))  
}
