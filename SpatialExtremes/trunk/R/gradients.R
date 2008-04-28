smithgrad <- function(par, data, distVec, loc.dsgn.mat, scale.dsgn.mat,
                      shape.dsgn.mat, fit.marge, std.err.type = "score"){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
    
  if (fit.marge){

     n.loccoeff <- ncol(loc.dsgn.mat)
     n.scalecoeff <- ncol(scale.dsgn.mat)
     n.shapecoeff <- ncol(shape.dsgn.mat)

    if (n.loccoeff == 1)
      loc.names <- "locCoeff"
    
    else
      loc.names <- paste("locCoeff", 1:n.loccoeff, sep="")
    
    if (n.scalecoeff == 1)
      scale.names <- "scaleCoeff"
    
    else
      scale.names <- paste("scaleCoeff", 1:n.scalecoeff, sep="")
    
    if (n.shapecoeff == 1)
      shape.names <- "shapeCoeff"
    
    else
      shape.names <- paste("shapeCoeff", 1:n.shapecoeff, sep="")
  
    param <- c("icov11", "icov12", "icov22", loc.names, scale.names, shape.names)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 5) == "scale")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    icov11 <- par["icov11"]
    icov12 <- par["icov12"]
    icov22 <- par["icov22"]
    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    icov11 <- par["icov11"]
    icov12 <- par["icov12"]
    icov22 <- par["icov22"]
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }
  
  grad <- .C("smithgrad", as.double(data), as.double(distVec), as.integer(n.site),
             as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
             as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
             as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
             as.double(shape.param), as.double(icov11), as.double(icov12),
             as.double(icov22), fit.marge, grad = double(n.obs * 3), PACKAGE = "SpatialExtremes")$grad

  grad <- matrix(grad, nrow = n.obs, ncol = 3)

  if (std.err.type == "score")
    jacobian <- var(grad) * n.obs

  if (std.err.type == "grad"){
    jacobian <- 0
    for (i in 1:n.obs){
      grad.vec <- matrix(grad[i,], ncol = 1)
      jacobian <- jacobian + grad.vec %*% t(grad.vec)
    }
  }

  return(jacobian)
}
