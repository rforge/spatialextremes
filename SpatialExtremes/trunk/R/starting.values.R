.start.smith <- function(data, coord, loc.model, scale.model, shape.model,
                         print.start.values = TRUE){

  cat("Computing appropriate starting values\n")
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    data[,i] <- gev2frech(data[,i], loc[i], scale[i], shape[i])
  }

  covs <- smithfull(data, coord, method = "Nelder")$param
  
  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0

  start <- as.list(covs)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))
  
  cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}


.start.schlather <- function(data, coord, loc.model, scale.model, shape.model,
                             print.start.values = TRUE){

  cat("Computing appropriate starting values\n")
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    data[,i] <- gev2frech(data[,i], loc[i], scale[i], shape[i])
  }

  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0

  cov.param <- schlatherfull(data, coord, method = "Nelder")$param
  
  start <- as.list(cov.param)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))
  
  cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}
