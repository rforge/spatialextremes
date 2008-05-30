TIC.maxstab <- function(object, ...){

  log.plik <- object$logLik
  ihessian <- object$ihessian
  jacobian <- object$jacobian

  if (is.null(ihessian) || is.null(jacobian))
    return(NA)

  penalty <- jacobian %*% ihessian
  
  tic <- -log.plik - sum(diag(penalty))
  return(tic)  
}
