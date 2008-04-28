##This functions fits the model with a penalized spline using radial
##basis functions for the GEV parameters. Thus, the coefficient of the
##polynomial models and the covariance matrix are estimated.
smithrb <- function(data, coord, knots, degree, penalty,
                    start, fit.marge = TRUE, ...,
                    warn.inf = TRUE, method = "BFGS",
                    std.err.type = "none", corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##locations is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pair <- n.site * (n.site - 1) / 2

  distVec <- distance(coord, vec = TRUE)
     
  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  loc.dsgn.mat <- modeldef.rb(coord, degree, knots, penalty)
  loc.pen.mat <- loc.dsgn.mat$pen.mat
  loc.dsgn.mat <- loc.dsgn.mat$dsgn.mat
  scale.dsgn.mat <- modeldef.rb(coord, degree, knots, penalty)
  scale.pen.mat <- scale.dsgn.mat$pen.mat
  scale.dsgn.mat <- scale.dsgn.mat$dsgn.mat
  shape.dsgn.mat <- modeldef.rb(coord, degree, knots, penalty)
  shape.pen.mat <- shape.dsgn.mat$pen.mat
  shape.dsgn.mat <- shape.dsgn.mat$dsgn.mat
  
  n.loccoeff <- ncol(loc.dsgn.mat)
  n.nonparloccoeff <- ncol(loc.pen.mat)
  n.scalecoeff <- ncol(scale.dsgn.mat)
  n.nonparscalecoeff <- ncol(scale.pen.mat)
  n.shapecoeff <- ncol(shape.dsgn.mat)
  n.nonparshapecoeff <- ncol(shape.pen.mat)

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

  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  body(nplk) <- parse(text = paste("-.C('smithrb', data, as.double(distVec), as.integer(n.site), as.integer(degree), as.integer(n.obs), loc.dsgn.mat, loc.pen.mat, as.integer(n.loccoeff), scale.dsgn.mat, scale.pen.mat, as.integer(n.scalecoeff), shape.dsgn.mat, shape.pen.mat, as.integer(n.shapecoeff), as.double(penalty),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                          "as.double(icov11), as.double(icov12), as.double(icov22), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk
  
  if (missing(start)) {

    locs <- scales <- rep(NA, n.site)
    shapes <- rep(0, n.site)

    dataFrech <- data
    
    for (i in 1:n.site){
      marg.param <- gevmle(data[,i])
      locs[i] <- marg.param["loc"]
      scales[i] <- marg.param["scale"]
      shapes[i] <- marg.param["shape"]
      dataFrech[,i] <- gev2frech(data[,i], locs[i], scales[i], shapes[i])
    }

    ##Compute the mean of each pairwise covariance matrix
    icovs <- NULL
    for (i in 1:(n.site-1))
      for (j in (i+1):n.site){
        icov.hat <- solve(cov(dataFrech[,c(i,j)]))
        icovs <- rbind(icovs, c(icov.hat[1,1], icov.hat[1,2],
                                icov.hat[2,2]))
      }

    icovs <- apply(icovs, 2, mean)
    meanCov <- solve(matrix(c(icovs[1], icovs[2], icovs[2],
                              icovs[3]), 2))
        
    start <- list(icov11 = meanCov[1,1], icov12 = meanCov[1,2],
                  icov22 = meanCov[2,2])

    locCoeff <- scaleCoeff <- shapeCoeff <- rep(0, n.loccoeff)

    locCoeff <- rbpspline(locs, coord, knots, degree, penalty)$beta
    scaleCoeff <- rbpspline(scales, coord, knots, degree, penalty)$beta
    shapeCoeff <- rbpspline(shapes, coord, knots, degree, penalty)$beta

    locCoeff[is.na(locCoeff)] <- 0
    scaleCoeff[is.na(scaleCoeff)] <- 0
    shapeCoeff[is.na(shapeCoeff)] <- 0

    names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL

    start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                          scaleCoeff = scaleCoeff,
                                          shapeCoeff = shapeCoeff))))
    
    start <- start[!(param %in% names(list(...)))]
  }

  if (!is.list(start)) 
    stop("`start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nplk)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nplk) <- c(f[m], f[-m])
  nllh <- function(p, ...) nplk(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  if (warn.inf && do.call("nllh", start.arg) == 1e36) 
    warning("negative log-likelihood is infinite at starting values")

  opt <- optim(start, nllh, hessian = hessian, ..., method = method)

  if ((opt$convergence != 0) || (opt$value == 1e36)) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
  }

  else opt$convergence <- "successful"

  param <- c(opt$par, unlist(fixed.param))
  
  iSigma <- matrix(c(param["icov11"], param["icov12"], param["icov12"],
                     param["icov22"]), 2, 2)

  extCoeff <- function(posVec)
    2 * pnorm(sqrt(posVec %*% iSigma %*% posVec) / 2)
  
  if (std.err.type == "observed"){
    
    tol <- .Machine$double.eps^0.5
    
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)){
      warning("observed information matrix is singular; passing std.err.type to ``none''")
      std.err.type <- "none"
      return
    }

    else{
      var.cov <- solve(var.cov, tol = tol)
      jacobian <- smithgrad(param, data, distVec, loc.dsgn.mat,
                            scale.dsgn.mat, shape.dsgn.mat,
                            fit.marge = fit.marge)
      
      var.cov <- var.cov[1:3,1:3] %*% jacobian %*% var.cov[1:3,1:3]
      
      std.err <- diag(var.cov)
      
      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        warning("Some (observed) standard errors are negative;\n passing them to NA")
        std.err[std.idx] <- NA
      }
      
      
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(var.cov) <- rownames(var.cov) <- 
        names(std.err) <- c("icov11", "icov12", "icov22")
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- NULL
  }
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data, est = "MLE",
                 logLik = -opt$value, opt.value = opt$value, model = "smith", coord = coord,
                 fit.marge = fit.marge, extCoeff = extCoeff, cov.mod = "Gaussian",
                 lik.fun = nllh)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}


##This functions fits the model with a polynomial response surface for
##the GEV parameters. Thus, the coefficient of the polynomial models
##and the covariance matrix are estimated.
smithlm <- function(data, coord, loc.form, scale.form, shape.form,
                    start, fit.marge = TRUE, ...,
                    warn.inf = TRUE, method = "BFGS",
                    std.err.type = "none", corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pair <- n.site * (n.site - 1) / 2

  distVec <- distance(coord, vec = TRUE)
     
  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  loc.dsgn.mat <- modeldef.lm(coord, loc.form)$dsgn.mat
  scale.dsgn.mat <- modeldef.lm(coord, scale.form)$dsgn.mat
  shape.dsgn.mat <- modeldef.lm(coord, shape.form)$dsgn.mat
  
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

  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  body(nplk) <- parse(text = paste("-.C('smithlm', data, as.double(distVec), as.integer(n.site), as.integer(n.obs), loc.dsgn.mat, as.integer(n.loccoeff), scale.dsgn.mat, as.integer(n.scalecoeff), shape.dsgn.mat, as.integer(n.shapecoeff),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                          "as.double(icov11), as.double(icov12), as.double(icov22), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk
  
  if (missing(start)) {

    dataFrech <- data
    if (fit.marge){
      locs <- scales <- rep(NA, n.site)
      shapes <- rep(0, n.site)
            
      for (i in 1:n.site){
        marg.param <- gevmle(data[,i])
        locs[i] <- marg.param["loc"]
        scales[i] <- marg.param["scale"]
        shapes[i] <- marg.param["shape"]
        dataFrech[,i] <- gev2frech(data[,i], locs[i], scales[i], shapes[i])
      }
    }

    ##Compute the mean of each pairwise covariance matrix
    meanCov <- 0
    for (i in 1:(n.site-1))
      for (j in (i+1):n.site)
        meanCov <- meanCov + cov(dataFrech[,c(i,j)]) / n.pair

    meanCov <- solve(meanCov)
    
    start <- list(icov11 = meanCov[1,1], icov12 = meanCov[1,2],
                  icov22 = meanCov[2,2])

    loc.form <- update(loc.form, locs ~ .)
    scale.form <- update(scale.form, scales ~ .)
    shape.form <- update(shape.form, shapes ~ .)
    
    locCoeff <- rep(0, n.loccoeff)
    scaleCoeff <- rep(0, n.scalecoeff)
    shapeCoeff <- rep(0, n.shapecoeff)
    locCoeff <- lm(loc.form, data = as.data.frame(cbind(locs, coord)))$coeff
    locCoeff[is.na(locCoeff)] <- 0
    scaleCoeff <- lm(scale.form, data = as.data.frame(cbind(scales, coord)))$coeff
    scaleCoeff[is.na(scaleCoeff)] <- 0
    shapeCoeff <- lm(shape.form, data = as.data.frame(cbind(shapes, coord)))$coeff
    shapeCoeff[is.na(shapeCoeff)] <- 0

    names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL

    start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                          scaleCoeff = scaleCoeff,
                                          shapeCoeff = shapeCoeff))))
    
    start <- start[!(param %in% names(list(...)))]
  }

  if (!is.list(start)) 
    stop("`start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nplk)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nplk) <- c(f[m], f[-m])
  nllh <- function(p, ...) nplk(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  if (warn.inf && do.call("nllh", start.arg) == 1e36) 
    warning("negative log-likelihood is infinite at starting values")

  opt <- optim(start, nllh, hessian = hessian, ..., method = method)

  if ((opt$convergence != 0) || (opt$value == 1e36)) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
  }

  else opt$convergence <- "successful"

  param <- c(opt$par, unlist(fixed.param))
  
  iSigma <- matrix(c(param["icov11"], param["icov12"], param["icov12"],
                     param["icov22"]), 2, 2)

  extCoeff <- function(posVec)
    2 * pnorm(sqrt(posVec %*% iSigma %*% posVec) / 2)
  
  if (std.err.type == "observed"){
    
    tol <- .Machine$double.eps^0.5
    
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)){
      warning("observed information matrix is singular; passing std.err.type to ``none''")
      std.err.type <- "none"
      return
    }

    else{
      var.cov <- solve(var.cov, tol = tol)
      jacobian <- smithgrad(param, data, distVec, loc.dsgn.mat,
                            scale.dsgn.mat, shape.dsgn.mat,
                            fit.marge = fit.marge)
      
      var.cov <- var.cov[1:3,1:3] %*% jacobian %*% var.cov[1:3,1:3]
      
      std.err <- diag(var.cov)
      
      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        warning("Some (observed) standard errors are negative;\n passing them to NA")
        std.err[std.idx] <- NA
      }
      
      
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(var.cov) <- rownames(var.cov) <- 
        names(std.err) <- c("icov11", "icov12", "icov22")
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- NULL
  }
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data, est = "MLE",
                 logLik = -opt$value, opt.value = opt$value, model = "smith", coord = coord,
                 fit.marge = fit.marge, extCoeff = extCoeff, cov.mod = "Gaussian",
                 loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                 lik.fun = nllh)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}

smithdsgnmat <- function(data, coord, loc.model, scale.model, shape.model,
                         start, fit.marge = TRUE, ..., warn.inf = TRUE,
                         method = "BFGS", std.err.type = "none", corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##locations is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pair <- n.site * (n.site - 1) / 2

  distVec <- distance(coord, vec = TRUE)
     
  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  ##Retrieve each element from loc.model, scale.model, shape.model
  loc.dsgn.mat <- loc.model$dsgn.mat
  scale.dsgn.mat <- scale.model$dsgn.mat
  shape.dsgn.mat <- shape.model$dsgn.mat

  loc.pen.mat <- loc.model$pen.mat
  scale.pen.mat <- scale.model$pen.mat
  shape.pen.mat <- shape.model$pen.mat

  loc.penalty <- loc.model$penalty
  scale.penalty <- scale.model$penalty
  shape.penalty <- shape.model$penalty

  ##Check if there some non penalizing terms
  if (is.null(loc.pen.mat))
    loc.pen.mat <- loc.penalty <- as.double(0.0)

  if (is.null(scale.pen.mat))
    scale.pen.mat <- scale.penalty <- as.double(0.0)

  if (is.null(shape.pen.mat))
    shape.pen.mat <- shape.penalty <- as.double(0.0)

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

  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  body(nplk) <- parse(text = paste("-.C('smithdsgnmat', data, as.double(distVec), as.integer(n.site), as.integer(n.obs), loc.dsgn.mat, loc.pen.mat, as.integer(n.loccoeff), as.double(loc.penalty), scale.dsgn.mat, scale.pen.mat, as.integer(n.scalecoeff), as.double(scale.penalty), shape.dsgn.mat, shape.pen.mat, as.integer(n.shapecoeff), as.double(shape.penalty),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                          "as.double(icov11), as.double(icov12), as.double(icov22), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk
  
  if (missing(start)) {

    start <- .start.smith(data, coord, loc.model, scale.model,
                          shape.model)
    
    start <- start[!(param %in% names(list(...)))]
  }

  if (!is.list(start)) 
    stop("`start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nplk)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nplk) <- c(f[m], f[-m])
  nllh <- function(p, ...) nplk(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  if (warn.inf && do.call("nllh", start.arg) == 1e36) 
    warning("negative log-likelihood is infinite at starting values")

  cat("Optimizing the negative log pairwise likelihood\n")
  opt <- optim(start, nllh, hessian = hessian, ..., method = method)
  cat("Optimization procedure is over\n")

  if ((opt$convergence != 0) || (opt$value == 1e36)) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
  }

  else opt$convergence <- "successful"

  param <- c(opt$par, unlist(fixed.param))
  
  iSigma <- matrix(c(param["icov11"], param["icov12"], param["icov12"],
                     param["icov22"]), 2, 2)

  extCoeff <- function(posVec)
    2 * pnorm(sqrt(posVec %*% iSigma %*% posVec) / 2)
  
  if (std.err.type == "observed"){
    
    tol <- .Machine$double.eps^0.5
    
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)){
      warning("observed information matrix is singular; passing std.err.type to ``none''")
      std.err.type <- "none"
      return
    }

    else{
      var.cov <- solve(var.cov, tol = tol)
      jacobian <- smithgrad(param, data, distVec, loc.dsgn.mat,
                            scale.dsgn.mat, shape.dsgn.mat,
                            fit.marge = fit.marge)
      
      var.cov <- var.cov[1:3,1:3] %*% jacobian %*% var.cov[1:3,1:3]
      
      std.err <- diag(var.cov)
      
      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        warning("Some (observed) standard errors are negative;\n passing them to NA")
        std.err[std.idx] <- NA
      }
      
      
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(var.cov) <- rownames(var.cov) <- 
        names(std.err) <- c("icov11", "icov12", "icov22")
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- NULL
  }
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data, est = "MLE",
                 logLik = -opt$value, opt.value = opt$value, model = "smith", coord = coord,
                 fit.marge = fit.marge, extCoeff = extCoeff, cov.mod = "Gaussian",
                 lik.fun = nllh, loc.model = loc.model, scale.model = scale.model,
                 shape.model = shape.model)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}


.was.modeldef.rb <- function(data, degree, knots, penalty){
  ##This function computes the design matrix X for using radial basis
  ##as well as the penalty matrix K

  X0 <- NULL
  for (i in 1:((degree - 1) / 2))
    X0 <- cbind(X0, data^i)
  
  X0 <- cbind(1, X0)
  
  ##Define the number of parametric parameters i.e. weights related to
  ##radial basis function are not taken into account
  n.par <- ncol(X0)
  
  X1 <- NULL

  if (is.null(dim(knots))){
    n.knots <- length(knots)

    for (k in 1:n.knots)
      X1 <- cbind(X1, as.matrix(dist(c(knots[k], data)))[-1,1])

  }

  else{
    n.knots <- nrow(knots)
    for (k in 1:n.knots)
      X1 <- cbind(X1, as.matrix(dist(rbind(knots[k,], data)))[-1,1])
  }   

  X1 <- X1^degree
  dsgn.mat <- cbind(X0, X1)

  pen.mat <- matrix(0, nrow = dim(dsgn.mat)[2],
                    ncol = dim(dsgn.mat)[2])
  pen.mat[-(1:n.par),-(1:n.par)] <- as.matrix(dist(knots, upper = TRUE, diag = TRUE))^degree

  dsgn.mat <- as.matrix(dsgn.mat)
  pen.mat <- as.matrix(pen.mat)

  init.fun <- function(y)
    rbpspline(y, data, knots, degree, penalty)$beta

  return(list(dsgn.mat = dsgn.mat, pen.mat = pen.mat, degree = degree,
              knots = knots, type = "rb", penalty.tot = penalty^degree,
              init.fun = init.fun, penalty = penalty))
}

##This functions fits the model with a polynomial response surface for
##the GEV parameters. Thus, the coefficient of the polynomial models
##and the covariance function parameters are estimated.
schlatherlm <- function(data, coord, cov.mod, loc.form, scale.form, shape.form,
                        start, fit.marge = TRUE, ...,
                        warn.inf = TRUE, method = "BFGS",
                        std.err.type = "none", corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pairs <- n.site * (n.site - 1) / 2
  
  dist <- distance(coord)
  
  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  if (!(cov.mod %in% c("whitmat","cauchy","powexp")))
    stop("``cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp'")

  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3

  loc.dsgn.mat <- modeldef.lm(coord, loc.form)$dsgn.mat
  scale.dsgn.mat <- modeldef.lm(coord, scale.form)$dsgn.mat
  shape.dsgn.mat <- modeldef.lm(coord, shape.form)$dsgn.mat
  
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
  
  param <- c("scale", "smooth", loc.names, scale.names, shape.names)

  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  body(nplk) <- parse(text = paste("-.C('schlatherlm', as.integer(cov.mod.num), data, as.double(dist), as.integer(dist.dim), as.integer(n.site), as.integer(n.obs), loc.dsgn.mat, as.integer(n.loccoeff), scale.dsgn.mat, as.integer(n.scalecoeff), shape.dsgn.mat, as.integer(n.shapecoeff),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                         "as.double(scale), as.double(smooth), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk
  
  if (missing(start)) {

    start <- list(scale = max(dist), smooth = .5)
    locs <- scales <- shapes <- rep(NA, n.site)
    
    for (i in 1:n.site){
      marg.param <- gevmle(data[,i])
      locs[i] <- marg.param["loc"]
      scales[i] <- marg.param["scale"]
      shapes[i] <- marg.param["shape"]
    }
        
    loc.form <- update(loc.form, locs ~ .)
    scale.form <- update(scale.form, scales ~ .)
    shape.form <- update(shape.form, shapes ~ .)
    
    locCoeff <- rep(0, n.loccoeff)
    scaleCoeff <- rep(0, n.scalecoeff)
    shapeCoeff <- rep(0, n.shapecoeff)
    locCoeff <- lm(loc.form, data = as.data.frame(cbind(locs, coord)))$coeff
    locCoeff[is.na(locCoeff)] <- 0
    scaleCoeff <- lm(scale.form, data = as.data.frame(cbind(scales, coord)))$coeff
    scaleCoeff[is.na(scaleCoeff)] <- 0
    shapeCoeff <- lm(shape.form, data = as.data.frame(cbind(shapes, coord)))$coeff
    shapeCoeff[is.na(shapeCoeff)] <- 0
    
    names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
    
    start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                          scaleCoeff = scaleCoeff,
                                          shapeCoeff = shapeCoeff))))
    start <- start[!(param %in% names(list(...)))]
  }
  

  if (!is.list(start)) 
    stop("`start' must be a named list")

  if (!length(start)) 
    stop("there are no parameters left to maximize over")

  nm <- names(start)
  l <- length(nm)
  f <- formals(nplk)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("`start' specifies unknown arguments")
  
  formals(nplk) <- c(f[m], f[-m])
  nllh <- function(p, ...) nplk(p, ...)
    
  if(l > 1)
    body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  if (warn.inf && do.call("nllh", start.arg) == 1e36) 
    warning("negative log-likelihood is infinite at starting values")

  opt <- optim(start, nllh, hessian = hessian, ..., method = method)

  if ((opt$convergence != 0) || (opt$value == 1e36)) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
  }

  else opt$convergence <- "successful"

  if (std.err.type == "observed"){
    
    tol <- .Machine$double.eps^0.5
    
    var.cov <- qr(opt$hessian, tol = tol)
    if(var.cov$rank != ncol(var.cov$qr)){
      warning("observed information matrix is singular; passing std.err.type to ``none''")
      std.err.type <- "none"
      return
    }

    else{
      var.cov <- solve(var.cov, tol = tol)
      
      std.err <- diag(var.cov)

      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        warning("Some (observed) standard errors are negative;\n passing them to NA")
        std.err[std.idx] <- NA
      }

      
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(var.cov) <- nm
      rownames(var.cov) <- nm
      names(std.err) <- nm
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- NULL
  }

  param <- c(opt$par, unlist(fixed.param))

  cov.fun <- covariance(scale = param["scale"], smooth = param["smooth"],
                        cov.mod = cov.mod)
  
  extCoeff <- function(h)
    1 + sqrt(1 - 1/2 * (cov.fun(h) + 1))

  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, param = param, cov.fun = cov.fun, fixed = unlist(fixed.param),
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, est = "MLE", data = data,
                 logLik = -opt$value, opt.value = opt$value, model = "schlather",
                 cov.mod = cov.mod, fit.marge = fit.marge, ext.coeff = extCoeff,
                 hessian = opt$hessian, lik.fun = nllh, loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form, coord = coord)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}
