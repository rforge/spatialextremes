gevmle <- function(x, method = "BFGS"){

  n <- length(x)
  
  nlgev <- function(param)
    ##param = c(loc, scale, shape)
    -.C("gevlik", as.double(x), as.integer(n),
        as.double(param[1]), as.double(param[2]),
        as.double(param[3]), dns = double(1),
        PACKAGE = "SpatialExtremes")$dns
  
  start <- c(loc = 0, scale = 0, shape = 0)
  start["scale"] <- sqrt(6 * var(x)) / pi
  start["loc"] <- mean(x) - 0.58 * start["scale"]
  
  opt <- optim(start, nlgev, method = method)
  param <- opt$par
  names(param) <- c("loc", "scale", "shape")

  return(param)
}


gpdmle <- function(x, threshold, method = "BFGS"){

  nlgpd <- function(param){
    ##param = c(scale, shape)
    -.C("gpdlik", as.double(exceed), as.integer(nat), as.double(threshold),
        as.double(param[1]), as.double(param[2]), dns = double(1),
        PACKAGE = "SpatialExtremes")$dns
  }

  high <- (x > threshold) & !is.na(x)
  exceed <- as.double(x[high])
  nat <- as.integer(exceed)

  start <- c(mean(exceed) - min(threshold), 0.0)
  
  opt <- optim(start, nlgpd, method = method)
  param <- opt$par
  names(param) <- c("scale", "shape")

  return(param)
}
