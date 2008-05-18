dpostmaxstab <- function(par, prior, cov.mod, data, coord,
                         loc.form, scale.form, shape.form, ...){
  ##Comptes the log-posterior density for a max-stable process
  ##data is a matrix with each column corresponds to one location
  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pairs <- n.site * (n.site - 1) / 2

  if (!(cov.mod %in% c("gauss", "whitmat","cauchy","powexp")))
    stop("``cov.mod'' must be one of 'gauss', 'whitmat', 'cauchy', 'powexp'")

  if (cov.mod == "gauss")
    distVec <- distance(coord, vec = TRUE)

  else
    dist <- distance(coord)
 
  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3

  ##With our notation, formula must be of the form y ~ xxxx
  loc.form <- update(loc.form, y ~ .)
  scale.form <- update(scale.form, y ~ .)
  shape.form <- update(shape.form, y ~ .)

  loc.model <- modeldef(coord, loc.form)
  scale.model <- modeldef(coord, scale.form)
  shape.model <- modeldef(coord, shape.form)

  loc.dsgn.mat <- loc.model$dsgn.mat
  scale.dsgn.mat <- scale.model$dsgn.mat
  shape.dsgn.mat <- shape.model$dsgn.mat

  loc.pen.mat <- loc.model$pen.mat
  scale.pen.mat <- scale.model$pen.mat
  shape.pen.mat <- shape.model$pen.mat

  loc.penalty <- loc.model$penalty.tot
  scale.penalty <- scale.model$penalty.tot
  shape.penalty <- shape.model$penalty.tot

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
  
  dprior <- do.call("dpriormaxstab", c(list(par = par), prior))

  if (cov.mod == "gauss"){
    smithlik <- function(par)
      .C("smithdsgnmat", data, as.double(distVec), as.integer(n.site),
         as.integer(n.obs), loc.dsgn.mat, loc.pen.mat, as.integer(n.loccoeff),
         as.integer(n.pparloc), as.double(loc.penalty), scale.dsgn.mat, scale.pen.mat,
         as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty),
         shape.dsgn.mat, shape.pen.mat, as.integer(n.shapecoeff), as.integer(n.pparshape),
         as.double(shape.penalty), as.double(par[3 + 1:n.loccoeff]),
         as.double(par[3 + n.loccoeff + 1:n.scalecoeff]),
         as.double(par[3 + n.loccoeff + n.scalecoeff + 1:n.shapecoeff]), as.double(par[1]),
         as.double(par[2]), as.double(par[3]), dns = double(1), PACKAGE = "SpatialExtremes")$dns

    postCont <- smithlik(par)
  }

  else{
    schlatherlik <- function(par)
      .C("schlatherdsgnmat", as.integer(cov.mod.num), data, as.double(dist),
         as.integer(dist.dim), as.integer(n.site), as.integer(n.obs), loc.dsgn.mat, loc.pen.mat,
         as.integer(n.loccoeff), as.integer(n.pparloc), as.double(loc.penalty),
         scale.dsgn.mat, scale.pen.mat, as.integer(n.scalecoeff), as.integer(n.pparscale),
         as.double(scale.penalty), shape.dsgn.mat, shape.pen.mat, as.integer(n.shapecoeff),
         as.integer(n.pparshape), as.double(shape.penalty), as.double(par[2 + 1:n.loccoeff]),
         as.double(par[2 + n.loccoeff + 1:n.scalecoeff]),
         as.double(par[2 + n.loccoeff + n.scalecoeff + 1:n.shapecoeff]), as.double(par[1]),
         as.double(par[2]), dns = double(1), PACKAGE = "SpatialExtremes")$dns

    postCont <- schlatherlik(par)
  }

  ans <- dprior + postCont
  
  return(ans)
}

dpriormaxstab <- function(par, mean, icov){
  par[1:2] <- log(par[1:2])
  as.double((par - mean) %*% icov %*% (par - mean))
}

gibbs <- function(n, init, prior, cov.mod, ...,
                  psd, thin, burn){

  np <- length(prior[[1]])
  param <- letters[1:np]
  dpst <- function(a) dpostmaxstab(a, prior, cov.mod, ...)

  mc <- .Call("gibbs", as.integer(n), as.integer(np),
              as.integer(thin), as.double(init), as.double(psd),
              quote(dpst(x)), env = new.env(), PACKAGE = "SpatialExtremes")
  
  naccept <- mc[[2]]
  nexternal <- mc[[3]]
  mc <- matrix(mc[[1]], ncol = np, byrow = TRUE)
  dimnames(mc) <- list(seq(0, n, thin), param)
  
  nexternal[np+1] <- sum(nexternal)/np
  naccept[np+1] <- sum(naccept)/np
  ar <- round(c(naccept/n, nexternal/n),2)
  ardn <- list(c("acc.rates","ext.rates"), c(param,"total"))
  ar <- matrix(ar, ncol = np+1, byrow = TRUE, dimnames = ardn)

  structure(mc, ar = ar)
}

prior <- function(mean, cov){

  if(mode(mean) != "numeric")
    stop("`mean' must be a numeric vector")
  
  if(!is.matrix(cov) || mode(cov) != "numeric")
    stop("`cov' must be a symmetric matrix")
  
  if(any(abs(cov - t(cov)) > .Machine$double.eps^0.5))
    warning("`cov' may not be symmetric")
  
  eg <- eigen(cov, symmetric = TRUE, only.values = TRUE)$values
  
  if(any(eg <= 0))
    warning("`cov' may not be positive definite")

  icov <- solve(cov)
  
  list(mean = mean, icov = icov)
}

