madogram <- function(data, coord, n.lag = 13,
                     gev.param = c(0, 1, 0),
                     which = c("mado", "ext"), ...){
  
  if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")
  
  n.site <- ncol(data)
  dist <- distance(coord)
  
  for (i in 1:n.site){
    param <- gevmle(data[,i])
    data[,i] <- gev2frech(data[,i], param[1], param[2],
                          param[3])
    data[,i] <- frech2gev(data[,i], gev.param[1],
                          gev.param[2], gev.param[3])
  }
  
  centers <- seq(0, max(dist), length = n.lag - 1)
  centers[1] <- centers[2] / 4
  bins <- c(0, (centers[-n.lag] + centers[1]) / 2, Inf)
  
  mado <- rep(0, length = n.lag)
  
  for (k in 1:n.lag){
    idx <- which((dist <= bins[k+1]) & (dist > bins[k]))

    if (length(idx)>0){
      site1 <- (idx-1) %/% n.site + 1
      site2 <- (idx-1) %% n.site + 1
      
      for (i in 1:length(idx))
        mado[k] <- mado[k] + sum(abs(data[,site1[i]] -
                                     data[,site2[i]]))
      
      mado[k] <- mado[k] / 2 / length(idx) / nrow(data)
    }
    
    else
      mado[k] <- NA
   
  }

  if (gev.param[3] == 0)
    ext.coeff <- exp(mado/gev.param[2])

  else
    ext.coeff <- gev2frech(gev.param[1] + mado / gamma(1 - gev.param[3]),
                           gev.param[1], gev.param[2], gev.param[3])

  if (length(which) == 2)
    par(mfrow=c(1,2))

  if (any(which == "mado"))
    plot(bins[-1], mado, ...)

  if (any(which == "ext"))
    plot(bins[-1], ext.coeff, ...)
  
  invisible(cbind(bins = bins[-1], madogram = mado, ext.coeff = ext.coeff))
}

extcoeff.emp <- function(data, coord, ..., prob = 0, plot = TRUE, lowess = TRUE){

   if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  z <- - 1 / log(prob)
  n.site <- ncol(data)

  dist <- distance(coord)

  ##First we need to transform data to unit Frechet using empirical CDF
  frech <- data
  for (i in 1:n.site){
    idx <- order(frech[,i])
    frech[idx,i] <- ppoints(frech[,i], a = 0)
  }

  frech <- - 1 / log(frech)
  x.bar <- colMeans(1/frech)

  lik.fun <- function(theta){
    frech.scaled <- apply(t(frech[,pair]) * x.bar[pair], 2, max)
    return(-sum(frech.scaled > z) * log(theta) + theta * sum(1/pmax(z, frech.scaled)))
  }

  ext.coeff <- rep(NA, n.site * (n.site - 1) / 2)

  k <- 1
  
  for (i in 1:(n.site - 1)){
    for (j in (i+1):n.site){
      pair <- c(i,j)
      ext.coeff[k] <- optimize(lik.fun, interval = c(1, 2))$minimum
      k <- k + 1
    }
  }

  if (plot){
    plot(dist, ext.coeff, ...)

    if (lowess)
      lines(lowess(dist, ext.coeff), ...)
  }
  
  invisible(cbind(dist = dist, ext.coeff = ext.coeff))
}
