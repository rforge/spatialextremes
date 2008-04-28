profile.maxstab <- function(fitted, param, range, n = 10,
                            plot = TRUE, ...){
  param.names <- names(fitted$fitted.values)
  n.param <- length(param.names)
  idx.prof <- which(param.names == param)
  fixed.param <- fitted$fixed
  fixed.names <- names(fixed.param)
  
  start <- as.list(fitted$fitted.values)
  start <- start[-idx.prof]

  fixed.values <- seq(range[1], range[2], length = n)

  llik <- rep(NA, n)
  par <- matrix(NA, ncol = n.param - 1, nrow = n)

  optfun <- fitted$lik.fun
  
  for (i in 1:n){
    fixed.val <- fixed.values[i]
    ##We need to modify the body of optfun for each fixed value
    if (length(fixed.names) > 0)
      body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-1),"]", collapse = ","),
                              ",", paste(fixed.names, "=", as.numeric(fixed.param),
                                         collapse = ","), ",", param, "=", fixed.val, ")"))

    else
      body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-1),"]", collapse = ","),
                              ",", param, "=", fixed.val, ")"))

  
    opt <- optim(start, optfun)
    llik[i] <- -opt$value
    par[i,] <- opt$par
  }

  ans <- cbind(fixed.values, llik, par)
  colnames(ans) <- c(param, "llik", param.names[-idx.prof])

  if (plot)
    plot(ans[,1], ans[,2], xlab = param, ylab = "log-likelihood", ...)
  
  return(ans)
  
}

profile2d.maxstab <- function(fitted, params, ranges, n = 10,
                              plot = TRUE, ...){
  param.names <- names(fitted$fitted.values)
  n.param <- length(param.names)
  idx.prof <- which(param.names %in% params)
  fixed.param <- fitted$fixed
  fixed.names <- names(fixed.param)
  
  start <- as.list(fitted$fitted.values)
  start <- start[-idx.prof]

  fixed.values1 <- seq(ranges[1,1], ranges[1,2], length = n)
  fixed.values2 <- seq(ranges[2,1], ranges[2,2], length = n)

  llik <- matrix(NA, n, n)

  optfun <- fitted$lik.fun
  
  for (i in 1:n){
    for (j in 1:n){
      fixed.val <- c(fixed.values1[i], fixed.values2[j])
      ##We need to modify the body of optfun for each fixed value
      if (length(fixed.names) > 0)
        body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-2),"]", collapse = ","),
                                ",", paste(fixed.names, "=", as.numeric(fixed.param),
                                           collapse = ","), ",",
                                paste(params, "=", fixed.val, collapse = ","), ")"))
      
      else
        body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-2),"]", collapse = ","),
                                ",", paste(params, "=", fixed.val, collapse = ","), ")"))
      
  
      llik[i,j] <- -optim(start, optfun)$value
    }
  }

  ans <- list(coord = cbind(fixed.values1, fixed.values2), llik = llik)
  colnames(ans$coord) <- c(params)

  if (plot)
    filled.contour(fixed.values1, fixed.values2, llik,
                   xlab = params[1], ylab = params[2], ...)

  return(ans)
  
}
