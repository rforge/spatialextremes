predict.maxstab <- function(object, newdata, ...){
  
  param <- object$param
  loc.form <- object$loc.form
  scale.form <- object$scale.form
  shape.form <- object$shape.form
  
  if (!missing(newdata)){
    data <- newdata
    
    if (is.null(dim(newdata)))
      data <- t(as.matrix(data))
  }

  else
    data <- object$coord

  loc.dsgnmat <- modeldef(data, loc.form)$dsgn.mat
  scale.dsgnmat <- modeldef(data, scale.form)$dsgn.mat
  shape.dsgnmat <- modeldef(data, shape.form)$dsgn.mat
                          
  idx.loc <- which(substr(names(param), 1, 8) == "locCoeff")
  idx.scale <- which(substr(names(param), 1, 10) == "scaleCoeff")
  idx.shape <- which(substr(names(param), 1, 10) == "shapeCoeff")

  loc.pred <- loc.dsgnmat %*% param[idx.loc]
  scale.pred <- scale.dsgnmat %*% param[idx.scale]
  shape.pred <- shape.dsgnmat %*% param[idx.shape]

  ans <- cbind(loc.pred, scale.pred, shape.pred)
  colnames(ans) <- c("loc", "scale", "shape")
  return(ans)
}

predict.pspline <- function(object, new.data, ...){

  if (missing(new.data))
    new.data <- object$x

  degree <- object$degree
  knots <- object$knots
  beta <- object$beta
  
  dsgn.mat <- rb(new.data, degree = degree, knots = knots,
                 penalty = NULL)$dsgn.mat

  y <- dsgn.mat %*% beta
  ans <- cbind(new.data, y)
  colnames(ans) <- c("x", "y.hat")
  
  return(ans)
}
