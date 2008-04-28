dpostmaxstab <- function(par, prior, model, cov.mod, data, coord,
                         loc.form, scale.form, shape.form, ...){
  ##Comptes the log-posterior density for a max-stable process
  ##data is a matrix with each column corresponds to one location
  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pairs <- n.site * (n.site - 1) / 2
  
  dist <- distance(coord)

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

  dprior <- do.call("dpriormaxstab", c(list(par = par), prior))

  if (model == "schlather"){
    schlather <- function(par)
      .C("schlatherlm", as.integer(cov.mod.num), as.double(data), as.double(dist),
         as.integer(dist.dim), as.integer(n.site), as.integer(n.obs),
         as.double(loc.dsgn.mat), as.integer(n.loccoeff), as.double(scale.dsgn.mat),
         as.integer(n.scalecoeff), as.double(shape.dsgn.mat), as.integer(n.shapecoeff),
         as.double(par[2 + 1:n.loccoeff]), as.double(par[2 + n.loccoeff + 1:n.scalecoeff]),
         as.double(par[2 + n.loccoeff + n.scalecoeff + 1:n.shapecoeff]), as.double(par[1]),
         as.double(par[2]), dns = double(1), PACKAGE = "SpatialExtremes")$dns

    postCont <- schlather(par)

    print(dprior)
    print(postCont)
    
    dprior + postCont
  }
}

dpriormaxstab <- function(par, mean, icov){
  par[1:2] <- log(par[1:2])
  as.numeric((par -mean) %*% icov %*% (par - mean))
}

gibbs <- function(n, init, prior, model, cov.mod, ...,
                  psd, thin, burn){

  np <- length(prior[[1]])
  param <- letters[1:np]
  dpst <- function(a) dpostmaxstab(a, prior, model, cov.mod,
                                   ...)

  mc <- .Call("gibbs", as.integer(n), as.integer(np),
              as.integer(thin), as.double(init), as.double(psd),
              quote(dpst(x)), new.env(), PACKAGE = "SpatialExtremes")
  
  naccept <- mc[[2]]
  nexternal <- mc[[3]]
  mc <- matrix(mc[[1]], ncol = np, byrow = TRUE)
  dimnames(mc) <- list(seq(0, n, thin), param)
  ##mc <- mc[-(1:burn), drop = FALSE]

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

