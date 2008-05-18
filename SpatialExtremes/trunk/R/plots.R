extcoeff <- function(fitted, n = 150, ...){

  model <- fitted$model
  extCoeff <- fitted$ext.coeff
  param <- fitted$param
  ##Define an appropriate range for the x-y axis
  if (model == "smith"){
    eps <- .Machine$double.eps
    det <- param["cov11"] * param["cov22"] - param["cov12"]^2
    denom <- param["cov22"] - 2 * param[2] + param["cov12"]
    x1 <- 2 * qnorm(1 - eps) / sqrt(param["cov22"] / det)
    x2 <- 2 * qnorm(1 - eps) / sqrt(param["cov11"] / det)
    x3 <-2 * qnorm(1 - eps) * sqrt(param["cov11"] * param["cov12"] / denom -
                               param["cov12"]^2 / denom)
      
    x.range <- 1.05 * c(-max(x1, x3), max(x1, x3))
    y.range <- 1.05 * c(-max(x2, x3), max(x2, x3))
  }

  if (model == "schlather"){
    fun <- function(h) abs(2 - extCoeff(h))
    init <- sqrt(sum(colMeans(fitted$coord)^2))
    opt1 <- optim(init, fun, method = "BFGS")$par

    y.range <- x.range <- 1.05 * c(-abs(opt1), abs(opt1))
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

  contour(xs, ys, extcoeff.hat, ...)
}
    

    
