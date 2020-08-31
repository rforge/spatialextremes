pmvnorm <- function(upper, mean = rep(0, dim), cov = diag(1, dim), n.mc = 1000){
    dim <- length(upper)

    if (!is.numeric(upper) || !is.vector(upper))
        stop("'upper' must be a numeric vector.")

    if (!is.numeric(cov) || !is.matrix(cov))
        stop("'cov' must be a numeric matrix!")

    if (any(dim(cov) != dim))
        stop('upper and cov does not match!')

    ## Convert the covariance matrix to a correlation matrix and
    ## update the upperBound accordingly
    upper <- (upper - mean) / sqrt(diag(cov))
    cov <- cov2cor(cov)

    tmp <-.C(C_pmvnorm2, as.integer(n.mc), as.integer(dim), as.double(cov),
             as.double(upper), est = double(1), err = double(1),
             PACKAGE = "SpatialExtremes")

    ans <- tmp$est
    attr(ans, "error") <- tmp$err
    return(ans)
}
