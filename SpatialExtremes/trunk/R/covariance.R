covariance <- function(fitted, scale, smooth, cov.mod = "whitmat",
                       plot = TRUE, dist, ...){

  if (!missing(fitted)){
    cov.mod <- fitted$cov.mod
    smooth <- fitted$param["smooth"]
    scale <- fitted$param["scale"]
  }

  if (cov.mod == "gauss")
    stop("``covariance'' is not implemented for the Smith's model")

  if (!(cov.mod %in% c("whitmat", "cauchy", "powexp")))
    stop("Invalid covariance model. ``cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp'")
  
  if (cov.mod == "whitmat"){
    if ((smooth <= 0) || (scale <= 0) || (smooth > 50))
      stop("invalid parameter for the whittle-matern covariance function")
    
    cov.fun <- function(dist) {
      idx <- which(dist == 0)
      
      ans <- 2^(1-smooth) / gamma(smooth) * (dist / scale)^smooth *
        besselK(dist / scale, smooth)
      ans[idx] <- 1
      return(ans)
    }
  }

  if (cov.mod == "cauchy"){
    if ((smooth <= 0) || (scale <= 0))
      stop("invalid parameter for the cauchy covariance function")
    
    cov.fun <- function(dist) (1 + (dist / scale)^2)^-smooth
  }

  if (cov.mod == "powexp"){
    if ((smooth < 0) || (smooth > 2) || (scale <= 0))
      stop("invalid parameter for the powered exponential covariance function")

    cov.fun <- function(dist) exp(-(dist / scale)^smooth)
  }

  if (plot){
    tmp.fun <- function(dist) cov.fun(dist) - 0.05
    xlimsup <- uniroot(tmp.fun, c(0, 10^16))$root
    plot(cov.fun, from = 0, to = xlimsup, ...)
  }

  if (!missing(dist)){
    cov.val <- cov.fun(dist)
    return(list(cov.fun = cov.fun, cov.val = cov.val))
  }

  else
    invisible(cov.fun)

  
}
