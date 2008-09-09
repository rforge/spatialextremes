extcoeff <- function(fitted, cov.mod, param,
                     n = 200, xlab, ylab, ...){

  if (all(class(fitted) != "maxstab"))
    stop("The 'extcoeff' function is only available for object of class 'maxstab'")
  
  if (missing(fitted) && (missing(cov.mod) || missing(param)))
    stop("You must specify either 'fitted' either 'cov.mod' AND 'param'")

  if (!missing(fitted)){
    if (ncol(fitted$coord) > 2)
      stop("It's not possible to use this function when the coordinate space has a dimension > 2")
  
    model <- fitted$model
    extCoeff <- fitted$ext.coeff
    param <- fitted$param
  }

  else{
    if (length(param) != 3)
      stop("You must specify all parameters in the Smith's or Schalther's model")
  
    if (cov.mod == "gauss"){
      model <- "Smith"
      names(param) <- c("cov11", "cov12", "cov22")
      Sigma <- matrix(c(param[1:2], param[2:3]), 2, 2)
      iSigma <- try(solve(Sigma), silent = TRUE)
      
      if (!is.matrix(iSigma) || (det(Sigma) <= 0) ||
          (param[1] <= 0))
        stop("You defined a non semi-definite positive matrix")
      
      extCoeff <- function(h)
        2 * pnorm(sqrt(h %*% iSigma %*% h) / 2)
    }
    
    else{
      model <- "Schlather"
      names(param) <- c("sill", "range", "smooth")
      extCoeff <- function(h)
        1 + sqrt(1 - .5 * (covariance(sill = param[1], range = param[2],
                                      smooth = param[3], cov.mod = cov.mod,
                                      plot = FALSE, dist = h)))
    }
  }
  
  ##Define an appropriate range for the x-y axis
  if (model == "Smith"){
    A <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                  param["cov22"]), 2, 2)
    
    eigen.values <- eigen(solve(A))
    eigen.vectors <- eigen.values$vectors
    eigen.values <- eigen.values$values

    r <- 25 / sqrt(2)
    axis1 <- eigen.vectors %*% c(r / sqrt(eigen.values[1]), 0)
    axis2 <- eigen.vectors %*% c(0, r / sqrt(eigen.values[2]))

    x.max <- max(abs(axis1[1]), abs(axis2[1]))
    y.max <- max(abs(axis1[2]), abs(axis2[2]))
    
    x.range <- 1.3 * c(-x.max, x.max)
    y.range <- 1.3 * c(-y.max, y.max)
    
  }

  if (model == "Schlather"){
    fun <- function(h) abs(1.838 - extCoeff(h))
    opt1 <- optim(param[2], fun, method = "BFGS")$par

    y.range <- x.range <- c(-abs(opt1), abs(opt1))
  }

  extcoeff.hat <- matrix(NA, nrow = n, ncol = n)

  xs <- seq(x.range[1], x.range[2], length = n)
  ys <- seq(y.range[1], y.range[2], length = n)

  if (model == "Smith"){
    for (i in 1:n)
      for (j in 1:n)
        extcoeff.hat[i,j] <- extCoeff(c(xs[i], ys[j]))
  }

  if (model == "Schlather"){
    for (i in 1:n)
      for (j in 1:n){
        h <- sqrt(xs[i]^2 + ys[j]^2)
        extcoeff.hat[i,j] <- extCoeff(h)
      }
  }

  if (missing(fitted))
    coord.names <- c("x", "y")

  else
    coord.names <- colnames(fitted$coord)

  if (missing(xlab))
    xlab <- coord.names[1]

  if (missing(ylab))
    ylab <- coord.names[2]
  
  contour(xs, ys, extcoeff.hat, xlab = xlab, ylab = ylab, ...)
}
    
map <- function(fitted, param = c("loc", "scale", "shape", "quant"),
                ..., ret.per = 100, ranges = apply(fitted$coord, 2, range),
                n = 80, col = terrain.colors(n), plot.contour = TRUE){

  if (ncol(fitted$coord) > 2)
    stop("It's not possible to use this function when the coordinate space has a dimension >= 2")
  
  x.range <- ranges[,1]
  y.range <- ranges[,2]

  ans <- matrix(NA, nrow = n, ncol = n)

  xs <- seq(x.range[1], x.range[2], length = n)
  ys <- seq(y.range[1], y.range[2], length = n)
  
  for (i in 1:n){
    new.data <- cbind(xs[i], ys)
    colnames(new.data) <- colnames(fitted$coord)
    param.hat <- predict(fitted, new.data)

    if (param == "loc")
      ans[,i] <- param.hat[,"loc"]
    
    if (param == "scale")
      ans[i,] <- param.hat[,"scale"]
  
    if (param == "shape")
      ans[i,] <- param.hat[,"shape"]
  
    if (param == "quant")
      ans[i,] <- .qgev(1 - 1/ret.per, param.hat[,"loc"],
                       param.hat[,"scale"], param.hat[,"shape"])
  }

  coord.names <- colnames(fitted$coord)
  xlab <- coord.names[1]
  ylab <- coord.names[2]
  
  image(xs, ys, ans, ..., col = col, xlab = xlab, ylab = ylab)

  if (plot.contour)
    contour(xs, ys, ans, add = TRUE)

  invisible(list(x = xs, y = ys, z = ans))
}
