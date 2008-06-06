madogram <- function(data, coord, n.lag = 100,
                     gev.param = c(0, 1, -1)){
  
  n.site <- ncol(data)
  dist <- distance(coord)

  for (i in 1:n.site){
    param <- gevmle(data[,i])
    data[,i] <- gev2frech(data[,i], param[1], param[2],
                          param[3])
    data[,i] <- frech2gev(data[,i], gev.param[1],
                          gev.param[2], gev.param[3])
  }

  lags <- seq(0, max(dist), length = n.lag)
  ans <- rep(0, length(lags))
  k <- 1

  for (lag in lags){
    idx <- which(dist <= lag)

    if (length(idx)>0){
      site1 <- (idx-1) %/% n.site + 1
      site2 <- (idx-1) %% n.site + 1
      
      for (i in 1:length(idx))
        ans[k] <- ans[k] + mean(abs(data[,site1[i]] -
                                    data[,site2[i]]))
      
      ans[k] <- ans[k] / 2 / length(idx)
    }

    else
      ans[k] <- 0

    k <- k+1
  }

  if (gev.param[3] == 0)
    ext.coeff <- exp(ans/gev.param[2])

  else
    ext.coeff <- gev2frech(gev.param[1] + ans / gamma(1 - gev.param[3]),
                           gev.param[1], gev.param[2], gev.param[3])

  ##par(mfrow=c(1,2))
  ##plot(lags, ans)
  points(lags, ext.coeff)
}
      
    
  
