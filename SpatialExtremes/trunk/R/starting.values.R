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
  ##names(covs) <- c("cov11", "cov12", "cov22")

  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0

  ##Then we transform observation to unit Frechet using the predict
  ##values for the GEV parameters
  ##loc <- loc.model$dsgn.mat %*% locCoeff
  ##scale <- scale.model$dsgn.mat %*% scaleCoeff
  ##shape <- shape.model$dsgn.mat %*% shapeCoeff
  ##
  ##for (i in 1:n.site)
  ##  data[,i] <- gev2frech(data[,i], loc[i], scale[i], shape[i])
  ##
  ##covs <- smithfull(data, coord)$param
  ##names(covs) <- c("cov11", "cov12", "cov22")

  ##start <- list(cov11 = covs["cov11"], cov12 = covs["cov12"],
  ##              cov22 = covs["cov22"])
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

  ##Then we transform observation to unit Frechet using the predict
  ##values for the GEV parameters
  ##loc <- loc.model$dsgn.mat %*% locCoeff
  ##scale <- scale.model$dsgn.mat %*% scaleCoeff
  ##shape <- shape.model$dsgn.mat %*% shapeCoeff
  ##
  ##for (i in 1:n.site)
  ##  data[,i] <- gev2frech(data[,i], loc[i], scale[i], shape[i])
  ##
  
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
