.define.parscale <- function(nllh, start, fixed.param, init.lik){
  n.par <- length(start)
  parscale <- rep(NA, n.par)
  
  for (i in 1:n.par){
    new.start <- start
    new.start[[i]] <- 1.1 * new.start[[i]]
    new.start.arg <- c(list(p = unlist(new.start)), fixed.param)
    parscale[i] <- do.call(nllh, new.start.arg) / init.lik
  }

  return(parscale)
  

}
