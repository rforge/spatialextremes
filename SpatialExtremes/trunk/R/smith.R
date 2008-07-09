##This file contains all the function to fit the max-stable
##characterisation of Smith


##This functions fits the model without any spatial structure for the
##GEV parameters. Thus, each GEV parameters are estimated at each
##location. However, if fit.marge = FALSE, observation are supposed to
##be unit Frechet and only the covariance matrix is estimated.
smithfull <- function(data, coord, start, fit.marge = FALSE,
                      ..., warn.inf = TRUE, method = "BFGS",
                      std.err.type = "none", control = list(),
                      corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pair <- n.site * (n.site - 1) / 2
  
  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  distVec <- distance(coord, vec = TRUE)
    
  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  if (fit.marge){

    loc.names <- paste("loc", 1:n.site, sep="")
    scale.names <- paste("scale", 1:n.site, sep="")
    shape.names <- paste("shape", 1:n.site, sep="")

    if (dist.dim == 2){
      param <- c("cov11", "cov12", "cov22", loc.names, scale.names, shape.names)
      body(nplk) <- parse(text = paste("-.C('smithfull', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs),",
                            paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                            paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                            paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                            "as.double(cov11), as.double(cov12), as.double(cov22), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
    }

    else{
      param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33", loc.names, scale.names, shape.names)
      body(nplk) <- parse(text = paste("-.C('smithfull3d', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs),",
                            paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                            paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                            paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                            "as.double(cov11), as.double(cov12), as.double(cov13), as.double(cov22), as.double(cov23), as.double(cov33), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
    }
  }
  
  else{

    if (dist.dim == 2){
      param <- c("cov11", "cov12", "cov22")
      body(nplk) <- parse(text = paste("-.C('smithfull', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs),",
                            paste("as.double(rep(1,", n.site, ")), "),
                            paste("as.double(rep(1,", n.site, ")), "),
                            paste("as.double(rep(1,", n.site, ")), "),
                            "as.double(cov11), as.double(cov12), as.double(cov22), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
    }

    else{
      param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")
      body(nplk) <- parse(text = paste("-.C('smithfull3d', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs),",
                            paste("as.double(rep(1,", n.site, ")), "),
                            paste("as.double(rep(1,", n.site, ")), "),
                            paste("as.double(rep(1,", n.site, ")), "),
                            "as.double(cov11), as.double(cov12), as.double(cov13), as.double(cov22), as.double(cov23), as.double(cov33), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
    }
  }
 
  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk

  if (missing(start)) {

    start <- NULL
    dataFrech <- data
    if (fit.marge){
      locs <- scales <- rep(NA, n.site)
      shapes <- rep(0, n.site)

      for (i in 1:n.site){
        marg.param <- gevmle(data[,i])
        locs[i] <- marg.param["loc"]
        scales[i] <- marg.param["scale"]
        shapes[i] <- marg.param["shape"]
      }

      start <- as.list(unlist(list(loc = locs, scale = scales, shape = shapes)))
    }

    exco.hat <- extcoeff.emp(data, coord, plot = FALSE)[,"ext.coeff"]
    sigma.init <- matrix(0, dist.dim, dist.dim)

    for (i in 1:n.pair)
      sigma.init <- solve(qnorm(exco.hat[i]/ 2)^2 *
                          (1 / distVec[i,]) %*%  t(1 / distVec[i,])) +
                            sigma.init

    print(sigma.init)
    sigma.init <- sigma.init / n.pair
    print(sigma.init)
    sigma.init <- solve(sigma.init)

    print(sigma.init)
    
    if (dist.dim == 2)
      start <- c(list(cov11 = sigma.init[1,1], cov12 = sigma.init[1,2],
                      cov22 = sigma.init[2,2]), start)

    else
      start <- c(list(cov11 = sigma.init[1,1], cov12 = sigma.init[1,2],
                      cov13 = sigma.init[1,3], cov22 = sigma.init[2,2],
                      cov23 = sigma.init[2,3], cov33 = sigma.init[3,3]), start)
    
    
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

  init.lik <- do.call("nllh", start.arg)
  if (warn.inf && (init.lik == 1.0e35)) 
    warning("negative log-likelihood is infinite at starting values")

  opt <- optim(start, nllh, hessian = hessian, ..., method = method,
               control = control)

  if ((opt$convergence != 0) || (opt$value == 1.0e35)) {
    warning("optimization may not have succeeded")

    if (opt$convergence == 1) 
      opt$convergence <- "iteration limit reached"
  }

  else opt$convergence <- "successful"

  if (opt$value == init.lik){
    warning("optimization stayed at the starting values. Consider tweaking the ndeps option.")
    opt$convergenc <- "Stayed at start. val."
  }

  param.names <- param
  param <- c(opt$par, unlist(fixed.param))
  
  if (std.err.type != "none"){
    
    var.cov <- try(solve(opt$hessian), silent = TRUE)

    if(!is.matrix(var.cov)){
      warning("observed information matrix is singular; passing std.err.type to ``none''")
      std.err.type <- "none"
      return
    }

    else{
      ihessian <- var.cov
      jacobian <- .smithgrad(param, data, distVec, as.double(0), as.double(0),
                             as.double(0), fit.marge = fit.marge, std.err.type =
                             std.err.type, fixed.param = names(fixed.param),
                             param.names = param.names)

      if(any(is.na(jacobian))){
        warning("observed information matrix is singular; passing std.err.type to ``none''")
        std.err.type <- "none"
        return
      }
      
      var.cov <- var.cov %*% jacobian %*% var.cov
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
    var.cov <- ihessian <- jacobian <- NULL
  }
    
  Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                    param["cov22"]), 2, 2)
  iSigma <- solve(Sigma)
  
  ext.coeff <- function(posVec)
    2 * pnorm(sqrt(posVec %*% iSigma %*% posVec) / 2)
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data, est = "MPLE",
                 logLik = -opt$value, opt.value = opt$value, model = "Smith",
                 fit.marge = fit.marge, ext.coeff = ext.coeff, cov.mod = "Gaussian",
                 lik.fun = nllh, coord = coord, ihessian = ihessian, jacobian = jacobian,
                 marg.cov = NULL)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}

##This functions fits the model from a generic R formula
##i.e. classical linear models as well as smoothing splines with
##radial basis functions may be used to model spatially the GEV
##parameters
smithform <- function(data, coord, loc.form, scale.form, shape.form,
                      start, fit.marge = TRUE, marg.cov = NULL, ...,
                      warn.inf = TRUE, method = "BFGS",
                      std.err.type = "none", control = list(),
                      corr = FALSE){
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

  ##With our notation, formula must be of the form y ~ xxxx
  loc.form <- update(loc.form, y ~ .)
  scale.form <- update(scale.form, y ~ .)
  shape.form <- update(shape.form, y ~ .)

  loc.model <- modeldef(cbind(coord, marg.cov), loc.form)
  scale.model <- modeldef(cbind(coord, marg.cov), scale.form)
  shape.model <- modeldef(cbind(coord, marg.cov), shape.form)

  loc.dsgn.mat <- loc.model$dsgn.mat
  scale.dsgn.mat <- scale.model$dsgn.mat
  shape.dsgn.mat <- shape.model$dsgn.mat

  loc.pen.mat <- loc.model$pen.mat
  scale.pen.mat <- scale.model$pen.mat
  shape.pen.mat <- shape.model$pen.mat

  loc.penalty <- loc.model$penalty.tot
  scale.penalty <- scale.model$penalty.tot
  shape.penalty <- shape.model$penalty.tot

  loc.type <- loc.model$type
  scale.type <- scale.model$type
  shape.type <- shape.model$type

  ##The total number of parameters to be estimated for each GEV
  ##parameter
  n.loccoeff <- ncol(loc.dsgn.mat)
  n.scalecoeff <- ncol(scale.dsgn.mat)
  n.shapecoeff <- ncol(shape.dsgn.mat)

  ##The number of ``purely parametric'' parameters to estimate i.e. we
  ##do not consider the weigths given to each basis function
  n.pparloc <- loc.model$n.ppar
  n.pparscale <- scale.model$n.ppar
  n.pparshape <- shape.model$n.ppar
  

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

  if (dist.dim == 2)
    param <- c("cov11", "cov12", "cov22", loc.names, scale.names, shape.names)

  else
    param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33", loc.names, scale.names, shape.names)
  
  param.names <- param

  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  if (dist.dim == 2)
    body(nplk) <- parse(text = paste("-.C('smithdsgnmat', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc), as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat), as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty), as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff), as.integer(n.pparshape), as.double(shape.penalty),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                          "as.double(cov11), as.double(cov12), as.double(cov22), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))

  else
    body(nplk) <- parse(text = paste("-.C('smithdsgnmat3d', as.double(data), as.double(distVec), as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc), as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat), as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty), as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff), as.integer(n.pparshape), as.double(shape.penalty),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                          "as.double(cov11), as.double(cov12), as.double(cov13), as.double(cov22), as.double(cov23), as.double(33), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
    
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

  init.lik <- do.call("nllh", start.arg)
  if (warn.inf && (init.lik == 1.0e35)) 
    warning("negative log-likelihood is infinite at starting values")

  if (is.null(control$parscale))
    control$parscale <- floor(as.numeric(start))

  control$parscale[control$parscale == 0] <- 1
  
  opt <- optim(start, nllh, hessian = hessian, ..., method = method,
               control = control)
  
  if ((opt$convergence != 0) || (opt$value == 1.0e35)) {
    warning("optimization may not have succeeded")

    if (opt$convergence != 0) 
      opt$convergence <- "iteration limit reached"
  }

  else opt$convergence <- "successful"

  if (opt$value == init.lik){
    warning("optimization stayed at the starting values. Consider tweaking the ndeps option.")
    opt$convergenc <- "Stayed at start. val."
  }

  param <- c(opt$par, unlist(fixed.param))
  
  if (std.err.type != "none"){
    
    var.cov <- try(solve(qr(opt$hessian)), silent = TRUE)
    if(!is.matrix(var.cov)){
      warning("observed information matrix is singular; passing std.err.type to ``none''")
      std.err.type <- "none"
      return
    }

    else{
      ihessian <- var.cov
      jacobian <- .smithgrad(param, data, distVec, loc.dsgn.mat,
                             scale.dsgn.mat, shape.dsgn.mat,
                             fit.marge = fit.marge, std.err.type =
                             std.err.type, fixed.param = names(fixed.param),
                             param.names = param.names)

      if(any(is.na(jacobian))){
        warning("observed information matrix is singular; passing std.err.type to ``none''")
        std.err.type <- "none"
        return
      }
      
      var.cov <- var.cov %*% jacobian %*% var.cov
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
        names(std.err) <- nm
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- ihessian <- jacobian <- NULL
  }

  if (dist.dim == 2)
    Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                      param["cov22"]), 2, 2)

  else
    Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov13"],
                      param["cov12"], param["cov22"], param["cov23"],
                      param["cov13"], param["cov23"], param["cov33"]),
                    3, 3)
  
  iSigma <- solve(Sigma)

  ext.coeff <- function(posVec)
    2 * pnorm(sqrt(posVec %*% iSigma %*% posVec) / 2)
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data, est = "MPLE",
                 logLik = -opt$value, opt.value = opt$value, model = "Smith", coord = coord,
                 fit.marge = fit.marge, ext.coeff = ext.coeff, cov.mod = "Gaussian",
                 loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                 lik.fun = nllh, loc.type = loc.type, scale.type = scale.type,
                 shape.type = shape.type, ihessian = ihessian, jacobian = jacobian,
                 marg.cov = marg.cov)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}
