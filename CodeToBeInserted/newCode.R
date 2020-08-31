rhitscenbrown <- function(n, coord, range, smooth, grid = FALSE){

    if (is.null(dim(coord))){
        dim <- 1
        n.site <- length(coord)
    } else {
        dim <- ncol(coord)
        n.site <- nrow(coord)
    }

    if (grid){
        n.eff.site <- n.site^dist.dim
        reg.grid <- .isregulargrid(coord[,1], coord[,2])
        steps <- reg.grid$steps
        reg.grid <- reg.grid$reg.grid
    }

    else
        n.eff.site <- n.site

    ans <- rep(-1e10, n * n.eff.site)
    ans2 <- rep(-1L, n * n.eff.site)


    ans <- .C(C_rhitscenbrown, as.double(coord), as.integer(n), as.integer(n.site), as.integer(dim),
              as.integer(grid), as.double(range), as.double(smooth), ans = ans, ans2 = ans2)$ans2

    if (grid){
        if (n == 1)
            ans <- matrix(ans, n.site, n.site)

        else
            ans <- array(ans, c(n.site, n.site, n))
    }

    else
        ans <- matrix(ans, nrow = n, ncol = n.site)

    return(ans)
}

pbrown <- function(data, coord, range, smooth, log = FALSE){
    n.site <- nrow(coord)
    n.obs <- nrow(data)
    dist.mat <- as.matrix(dist(coord))
    vario.mat <- (dist.mat / range)^smooth
    ## Computation of the V part: (see the formula of Huser and Davison in Biometrika for instance)
    V.contrib <- rep(0, n.obs)
    for (j in 1:n.site){
        dummy <- matrix(vario.mat[j,-j], n.site-1, n.site-1)
        Sigma <- (dummy + t(dummy) - vario.mat[-j,-j]) / (2 * sqrt(dummy * t(dummy)))

        for (i in 1:n.obs){
            mu <- sqrt(0.5 * vario.mat[j,-j]) - log(data[i,j] / data[i,-j]) / sqrt(2 * vario.mat[j,-j])
            V.contrib[i] <- V.contrib[i] + pmvnorm(upper = mu, cov = Sigma) / data[i,j]
        }
    }

    ans <- -V.contrib
    if (!log)
        ans <- exp(ans)

    return(ans)
}

efffitbrown <- function(data, coord, init, n.mc = 10^4, cutoff = 0.9, n.hit = 10, weighted = FALSE, n.mc2 = 100){

    ##data is a matrix with each column corresponds to one location
    ##locations is a matrix giving the coordinates (1 row = 1 station)
    n.site <- ncol(data)
    n.obs <- nrow(data)
    dist.dim <- ncol(coord)

    dist.mat <- as.matrix(dist(coord))

    nllik <- function(param){
        range <- param[1]
        smooth <- param[2]

        if (range <= 0)
            return((1 - range) * 1e6)

        if (smooth <= 0)
            return((1 - smooth) * 1e6)

        if (smooth > 2)
            return((smooth - 1) * 1e6)

        vario.mat <- (dist.mat / range)^smooth
        vario.vec <- (sqrt(rowSums(coord^2)) / range)^smooth

        ## log-density is of the form -V(z) + log(sum_tau prod_{j=1}^|tau| lambda_j * proba_j)

        ## Computation of the exp(-V) part: (see the formula of Huser and Davison in Biometrika for instance)
        log.proba <- pbrown(data, coord, range, smooth, log = TRUE)

        ## Now use the approximated part of the sum over all possible partition

        ## We need first to simulate a bunch of hitting scenario and
        ## estimate its empirical distribution (parameters being fixed)
        ## hit.scen <- rhitscenbrown(n.mc, coord, range, smooth)
        ## emp.freq <- sort(table(apply(hit.scen, 1, paste, collapse = "-")) / n.mc, decreasing = TRUE)
        ## idx <- which(cumsum(emp.freq) < cutoff)
        ##
        ## if (length(idx) == 0)
        ##     emp.freq <- emp.freq
        ## else if (length(idx) > n.hit)
        ##     emp.freq <- emp.freq[1:n.hit]
        ## else
        ##     emp.freq <- emp.freq[cumsum(emp.freq) < cutoff]
        ##
        ## ## Then compute the contribution
        ##
        ## ## Computation of fixed quantities
        ##
        ## Sigma <- -vario.mat## <<-- will be the covariance matrix
        ## for (i in 1:n.site)
        ##     for (j in 1:n.site)
        ##         Sigma[i,j] <- Sigma[i,j] + vario.vec[i] + vario.vec[j]
        ##
        ## iSigma <- solve(Sigma)
        ## one <- rep(1, n.site)
        ## mahal1 <- as.numeric(t(one) %*% iSigma %*% one)
        ## mahalVario <- as.numeric(t(vario.vec) %*% iSigma %*% vario.vec)
        ## tmp <- iSigma %*% one
        ## Q.mat <- iSigma - tmp %*% t(tmp) / mahal1
        ## tmp2 <- as.numeric(t(tmp) %*% vario.vec - 1)
        ## L.vec <- t(tmp2 / mahal1 * one - vario.vec) %*%
        ##     iSigma
        ##
        ## idx.tau <- 1
        ## overall.tau.contrib <- rep(0, n.obs)
        ## for (tau in names(emp.freq)){
        ##     tau <- as.numeric(unlist(strsplit(tau, split = "-")))
        ##     tau.size <- max(tau)
        ##
        ##     tau.contrib <- rep(1, n.obs)
        ##     for (label in 1:tau.size){
        ##         ## Computation of the intensity part, i.e., lambda_j
        ##         idx <- which(tau == label)
        ##         vario.vec.label <- vario.vec[idx]
        ##         Sigma.label <- Sigma[idx,idx,drop=FALSE]
        ##         iSigma.label <- solve(Sigma.label)
        ##         one.label <- rep(1, length(idx))
        ##
        ##         mahal1.label <- as.numeric(t(one.label) %*% iSigma.label %*% one.label)
        ##         mahalVario.label <- as.numeric(t(vario.vec.label) %*% iSigma.label %*% vario.vec.label)
        ##         tmp.label <- iSigma.label %*% one.label
        ##         Q.mat.label <- iSigma.label - tmp.label %*% t(tmp.label) / mahal1.label
        ##         tmp2.label <- as.numeric(t(tmp.label) %*% vario.vec.label - 1)
        ##         L.vec.label <- t(tmp2.label / mahal1.label * one.label - vario.vec.label) %*%
        ##             iSigma.label
        ##         C.label <- (2 * pi)^(0.5 * (1 - length(idx))) / sqrt(det(Sigma.label) * mahal1.label) *
        ##             exp(0.5 * tmp2.label^2 / mahal1.label - 0.5 * mahalVario.label)
        ##
        ##         if (length(idx) < n.site)
        ##             Sigma.cond <- solve(Q.mat[-idx,-idx])##The covariance mat for the proba part
        ##
        ##         for (i in 1:n.obs){
        ##             obs.label <- data[i,idx]
        ##             lambda.label <- C.label * exp(-0.5 * as.numeric(t(log(obs.label)) %*% Q.mat.label %*% log(obs.label)) +
        ##                                           L.vec.label %*% log(obs.label)) / prod(obs.label)
        ##
        ##             ## Computation of the proba part, i.e., proba_j ## Watch out overwrite existing quantities
        ##             if (length(idx) < n.site){
        ##                 mu.cond <- as.numeric((L.vec[-idx] - log(obs.label) %*% Q.mat[idx,-idx]) %*% Sigma.cond)
        ##                 proba.label <- pmvnorm(upper = log(data[i,-idx]), mean = mu.cond, cov = Sigma.cond, n.mc = n.mc2)
        ##             } else
        ##                 proba.label <- 1
        ##
        ##             tau.contrib[i] <- tau.contrib[i] * lambda.label * proba.label
        ##         }
        ##     }
        ##
        ##     if (weighted)
        ##         overall.tau.contrib <- overall.tau.contrib + tau.contrib * emp.freq[idx.tau]
        ##     else
        ##         overall.tau.contrib <- overall.tau.contrib + tau.contrib
        ##
        ##     idx.tau <- idx.tau + 1
        ##
        ## }
        ##
        ## ans <- -sum(log.proba + log(overall.tau.contrib))

        if (!is.finite(ans))
            browser()

        return(ans)
    }

    if (missing(init))
        init <- fitcovariance(data, coord, "brown")$param

    opt <- nlm(nllik, init, hessian = TRUE, print.level = 2, ndigit = 3)
    ##opt <- optim(init, nllik, hessian = TRUE, control = list(trace = 5), method = "BFGS")
    opt$counts <- opt$iter
    names(opt$counts) <- "function"
    opt$value <- opt$minimum
    opt$par <- opt$estimate
    var.cov <- solve(opt$hessian)
    corr.mat <- cov2cor(var.cov)
    std.err <- sqrt(diag(var.cov))
    names(opt$par) <- c("range", "smooth")
    param <- opt$par

    ext.coeff <- function(h)
        2 * pnorm((h / param["range"])^(0.5 * param["smooth"]) / sqrt(2))

    conc.prob <- function(h){
        n.sim <- 1000
        semivario <- matrix((h / param["range"])^param["smooth"], 2 * n.sim, length(h), byrow = TRUE)
        eps <- rnorm(n.sim)
        u1 <- pnorm(eps)
        eps <- c(eps, -eps)## antithetic
        u1 <- c(u1, 1 - u1)

        colMeans(1 / (u1 + exp(semivario - sqrt(2 * semivario) * eps) *
                          pnorm(sqrt(2 * semivario) - eps)))

    }

    fitted <- list(fitted.values = opt$par, std.err = std.err,
                   var.cov = var.cov, param = opt$par, cov.fun = NA, fixed = NULL,
                   deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                   counts = opt$counts, message = opt$message, est = "MLE", data = data,
                   logLik = -opt$value, opt.value = opt$value, model = "Brown-Resnick",
                   cov.mod = "brown", fit.marge = FALSE, ext.coeff = ext.coeff, iso = TRUE,
                   hessian = opt$hessian, lik.fun = nllik, coord = coord, ihessian = var.cov,
                   marg.cov = NULL, nllh = nllik, weighted = weighted, conc.prob = conc.prob)

    class(fitted) <- c(fitted$model, "maxstab")
    return(fitted)
}
