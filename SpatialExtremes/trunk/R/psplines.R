rbpspline <- function(y, x, knots, degree, penalty, ...){

  if ((degree %% 2) == 0)
    stop("``degree'' must be an odd number")

  if (missing(penalty)){
    cv.fit <- cv(y, x, knots, degree, ...)
    penalty <- cv.fit$minimum
    cv.hat <- cv.fit$objective
  }

  else
    cv.hat <- NA

  dsgn.mat <- rb(x, degree = degree, knots = knots,
                 penalty = penalty)
  pen.mat <- dsgn.mat$pen.mat
  dsgn.mat <- dsgn.mat$dsgn.mat

  ##To compute the beta vector, we use the Demmler-Reinch
  ##orthogonalization
  chol.mat <- chol(crossprod(dsgn.mat))
  chol.mat.inv <- solve(chol.mat)
  svd.mat <- svd(t(chol.mat.inv) %*% pen.mat %*% chol.mat.inv)

  A <- dsgn.mat %*% chol.mat.inv %*% svd.mat$u
  b <- t(A) %*% y

  fitted.values <- A %*% (b / (1 + penalty^degree * svd.mat$d))
  smooth.mat <- A %*% diag(1 / (1 + penalty^degree * svd.mat$d)) %*% t(A)
  beta <- chol.mat.inv %*% svd.mat$u %*% (b / (1 + penalty^degree * svd.mat$d))
  
  df <- sum(diag(smooth.mat))
  res.df <- length(y) - 2 * sum(diag(smooth.mat)) +
    sum(diag(smooth.mat %*% t(smooth.mat)))
  rank <- sum(eigen(smooth.mat, only.values = TRUE)$values != 0)
  
  
  fitted <- list(x = x, y = y, penalty = penalty, knots = knots,
                 degree = degree, fitted.values = fitted.values,
                 dsgn.mat = dsgn.mat, pen.mat = pen.mat, cv = cv.hat,
                 smooth.mat = smooth.mat, df = df, res.df = res.df,
                 rank = rank, beta = beta, call = match.call())

  class(fitted) <- "pspline"
  return(fitted)
}

cv <- function(y, x, knots, degree, pen.range = c(0, 1000),
               ...){

  if (any(pen.range < 0))
    stop("penalty term cannot be negative. Set an appropriate ``pen.range''")

  dsgn.mat <- rb(x, degree = degree, knots = knots, penalty = NULL)
  pen.mat <- dsgn.mat$pen.mat
  dsgn.mat <- dsgn.mat$dsgn.mat

  ##To compute the beta vector, we use the Demmler-Reinch
  ##orthogonalization
  chol.mat <- chol(crossprod(dsgn.mat))
  chol.mat.inv <- solve(chol.mat)
  svd.mat <- svd(t(chol.mat.inv) %*% pen.mat %*% chol.mat.inv)

  A <- dsgn.mat %*% chol.mat.inv %*% svd.mat$u
  b <- t(A) %*% y
  
  obj.fun <- function(penalty){
    fitted.values <- A %*% (b / (1 + penalty^degree * svd.mat$d))
    smooth.mat <- A %*% diag(1 / (1 + penalty^degree * svd.mat$d)) %*% t(A)
    
    cv <- sum(((y - fitted.values) / (1 - diag(smooth.mat)))^2)
    return(cv)
  }

  opt <- optimize(obj.fun, interval = pen.range, ...)

  if (opt$value == 1000)
    warning("the smoothing parameter estimate is equal to the upper bound.
You should change pen.range")
  return(opt)    

}

rb <- function(..., knots, degree, penalty){

  ##vars <- as.list(substitute(list(...)))[-1]

  ##data <- NULL
  ##for (i in 1:length(vars))
  ##  data <- cbind(data, eval(vars[[i]]))

  data <- cbind(...)
    
  if ((degree %% 2) == 0)
    stop("``degree'' must be an odd number")
  
  if (missing(penalty))
    penalty <- NULL

  X0 <- NULL
  for (i in 1:((degree - 1) / 2))
    X0 <- cbind(X0, data^i)
  
  X0 <- cbind(1, X0)
  
  ##Define the number of ``purely parametric'' parameters to be
  ##estimated i.e. weights related to radial basis function are not
  ##taken into account
  n.ppar <- ncol(X0)
  
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
  pen.mat[-(1:n.ppar),-(1:n.ppar)] <- as.matrix(dist(knots, upper = TRUE, diag = TRUE))^degree

  dsgn.mat <- as.matrix(dsgn.mat)
  pen.mat <- as.matrix(pen.mat)

  ans <- list(dsgn.mat = dsgn.mat, pen.mat = pen.mat, degree = degree,
              penalty = penalty, knots = knots, data = data, n.ppar = n.ppar,
              call = match.call())

  return(ans)
}
