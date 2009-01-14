.start.smith <- function(data, coord, covariables, loc.form, scale.form, shape.form,
                         print.start.values = TRUE, method = "Nelder",
                         iso = TRUE, ...){

  n.site <- ncol(data)
  
  if (ncol(coord) == 2)
    param <- c("cov11", "cov12", "cov22")

  else
    param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")

  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, marge = "emp", iso = iso),
              fixed.param)
    covs <- do.call("fitcovmat", args)$param
  }

  else
    covs <- fitcovmat(data, coord, marge = "emp", iso = iso)$param

  if (iso){
    covs <- covs[1]
    names(covs) <- "cov"
  }    

  spatgev <- fitspatgev(data, as.matrix(covariables), loc.form,
                          scale.form, shape.form, std.err.type = "none")

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  covs <- smithfull(frech, coord, fit.marge = FALSE, start = as.list(covs),
                    iso = iso, warn = FALSE, method = method, std.err.type = "none")$param
  
  start <- c(as.list(covs), as.list(spatgev$param))

  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(unlist(start))
  }

  return(start)
}


.start.schlather <- function(data, coord, covariables, cov.mod, loc.form,
                             scale.form, shape.form, print.start.values = TRUE,
                             method = "Nelder",
                             ...){

  n.site <- ncol(data)
  param <- c("sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  spatgev <- fitspatgev(data, as.matrix(covariables), loc.form,
                        scale.form, shape.form, std.err.type = "none")

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  cov.param <- schlatherfull(frech, coord, cov.mod = cov.mod, fit.marge = FALSE,
                             warn = FALSE, start = as.list(cov.param), method = method,
                             std.err.type = "none")$param
  
  start <- c(as.list(cov.param), as.list(spatgev$param))

  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}


.start.schlatherind <- function(data, coord, covariables, cov.mod, loc.form,
                                scale.form, shape.form, print.start.values = TRUE,
                                method = "Nelder", ...){

  n.site <- ncol(data)
  param <- c("alpha", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp", sill = 1), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  spatgev <- fitspatgev(data, as.matrix(covariables), loc.form,
                          scale.form, shape.form, std.err.type = "none")

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  cov.param <- schlatherindfull(frech, coord, cov.mod = cov.mod, fit.marge = FALSE,
                                start = c(list(alpha = 0.25), as.list(cov.param)),
                                warn = FALSE, method = method, std.err.type = "none")$param
  
  start <- c(as.list(cov.param), as.list(spatgev$param))
  
  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}

.start.geomgauss <- function(data, coord, covariables, cov.mod, loc.form,
                             scale.form, shape.form, print.start.values = TRUE,
                             method = "Nelder", ...){

  n.site <- ncol(data)
  param <- c("sigma2", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  spatgev <- fitspatgev(data, as.matrix(covariables), loc.form,
                          scale.form, shape.form, std.err.type = "none")

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  cov.param <- geomgaussfull(frech, coord, cov.mod = cov.mod, fit.marge = FALSE,
                             start = c(list(sigma2 = 1), as.list(cov.param)),
                             warn = FALSE, method = method, std.err.type = "none")$param
  
  start <- c(as.list(cov.param), as.list(spatgev$param))
  
  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}

.start.nsgeomgauss <- function(data, coord, covariables, cov.mod, loc.form,
                               scale.form, shape.form, sigma2.form,
                               print.start.values = TRUE, method = "Nelder",
                               ...){

  n.site <- ncol(data)
  sigma2.terms <- terms(sigma2.form)

  if (length(attributes(sigma2.terms)$factor) == 0)
    n.sigma2coeff <- attributes(sigma2.terms)$intercept

  else
    n.sigma2coeff <- attributes(sigma2.terms)$intercept +
      ncol(attributes(sigma2.terms)$factors)

  if (n.sigma2coeff == 1)
    sigma2.names <- "sigma2"

  else
    sigma2.names <- paste("sigma2Coeff", 1:n.sigma2coeff, sep="")

  param <- c(sigma2.names, "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  spatgev <- fitspatgev(data, as.matrix(covariables), loc.form,
                          scale.form, shape.form, std.err.type = "none")

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  sigma2param <- rep(0, length(sigma2.names))
  names(sigma2param) <- sigma2.names
  
  cov.param <- nsgeomgaussfull(frech, coord, cov.mod = cov.mod, fit.marge = FALSE,
                               sigma2.form = sigma2.form, warn = FALSE, method = method,
                               start = c(as.list(sigma2param), as.list(cov.param)),
                               std.err.type = "none")$param
  
  start <- c(as.list(cov.param), as.list(spatgev$param))

  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}

