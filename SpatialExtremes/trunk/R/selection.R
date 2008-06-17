TIC.maxstab <- function(object, ...){

  all.objects <- c(list(object), list(...))
  tic <- NULL
  for (i in 1:length(all.objects)){

    object <- all.objects[[i]]
    log.plik <- object$logLik
    ihessian <- object$ihessian
    jacobian <- object$jacobian
    
    if (is.null(ihessian) || is.null(jacobian))
      tic <- c(tic, NA)
    
    else{
      penalty <- jacobian %*% ihessian
      
      tic <- c(tic, -log.plik - sum(diag(penalty)))
    }
  }

  names(tic) <- as.character(sys.call())[-1]
  tic <- tic[order(tic)]
  
  return(tic)  
}
