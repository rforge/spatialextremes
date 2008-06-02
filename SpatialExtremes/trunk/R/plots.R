extcoeff <- function(fitted, n = 150, ...){

  model <- fitted$model
  extCoeff <- fitted$ext.coeff
  param <- fitted$param
  ##Define an appropriate range for the x-y axis
  if (model == "smith"){
    eps <- .Machine$double.eps
    det <- param["cov11"] * param["cov22"] - param["cov12"]^2
    denom1 <- param["cov22"] - 2 * param["cov12"] + param["cov11"]
    denom2 <- param["cov22"] + 2 * param["cov12"] + param["cov11"]
    x.lim <- 2 * qnorm(1 - eps) / sqrt(param["cov22"] / det)
    y.lim <- 2 * qnorm(1 - eps) / sqrt(param["cov11"] / det)
    x.diag1 <- 2 * qnorm(1 - eps) * sqrt(det / denom1)
    x.diag2 <- 2 * qnorm(1 - eps) * sqrt(det / denom2)
      
    x.range <- 1.2 * c(-max(x.lim, x.diag1, x.diag2), max(x.lim, x.diag1, x.diag2))
    y.range <- 1.2 * c(-max(y.lim, x.diag1, x.diag2), max(y.lim, x.diag1, x.diag2))
  }

  if (model == "schlather"){
    fun <- function(h) abs(2 - extCoeff(h))
    init <- sqrt(sum(colMeans(fitted$coord)^2))
    opt1 <- optim(init, fun, method = "BFGS")$par

    y.range <- x.range <- 1.2 * c(-abs(opt1), abs(opt1))
  }

  extcoeff.hat <- matrix(NA, nrow = n, ncol = n)

  xs <- seq(x.range[1], x.range[2], length = n)
  ys <- seq(y.range[1], y.range[2], length = n)

  if (model == "smith"){
    for (i in 1:n)
      for (j in 1:n)
        extcoeff.hat[i,j] <- extCoeff(c(xs[i], ys[j]))
  }

  if (model == "schlather"){
    for (i in 1:n)
      for (j in 1:n){
        h <- sqrt(xs[i]^2 + ys[j]^2)
        extcoeff.hat[i,j] <- extCoeff(h)
      }
  }

  coord.names <- colnames(fitted$coord)
  xlab <- coord.names[1]
  ylab <- coord.names[2]
  contour(xs, ys, extcoeff.hat, xlab = xlab, ylab = ylab, ...)
}
    
map <- function(fitted, param = c("loc", "scale", "shape", "quant"),
                ..., ret.per = 100, ranges = apply(fitted$coord, 2, range),
                n = 80, col = terrain.colors(n), plot.contour = TRUE){

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
