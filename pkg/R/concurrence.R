concprob <- function(data, coord, fitted, n.bins, add = FALSE, xlim = c(0, max(dist)),
                     ylim = c(min(0, concProb), max(1, concProb)), col = 1:2,
                     which = "kendall", block.size = floor(sqrt(nrow(data))),
                     plot = TRUE, ...){

    if (missing(fitted) && missing(data) && missing(coord))
        stop("You must either specify a fitted model OR 'data' and 'coord'")

    fit.curves <- FALSE
    if (!missing(fitted)){
        data <- fitted$data
        coord <- fitted$coord

        if ((fitted$model == "Smith") && !fitted$iso)
            warning("The concurrence probability function is valid only for isotropic models. The fitted curves won't be plotted.")

        else {
            fit.curves <- TRUE
            conc.prob.fit <- fitted$conc.prob
        }
    }

    if (is.null(dim(coord))){
        if (length(coord) != ncol(data))
            stop("'data' and 'coord' don't match")
    }

    else if (nrow(coord) != ncol(data))
        stop("'data' and 'coord' don't match")


    if (any(!(which %in% c("kendall", "emp", "boot"))))
        stop("'which' must be either 'kendall', 'emp' or 'boot'")
    
    n.obs <- nrow(data)
    n.site <- ncol(data)
    n.pairs <- n.site * (n.site - 1) / 2
    n.block <- floor(n.obs / block.size)
    dist <- distance(coord)

    ## Compute p(x1, x2) for each pair of station
    if (which == "kendall")
        concProb <- .C("concProbKendall", as.double(data), as.integer(n.site),
                       as.integer(n.obs), concProb = double(n.pairs), NAOK = TRUE,
                       PACKAGE = "SpatialExtremes")$concProb
    
    else if (which == "emp")
        concProb <- .C("empiricalConcProb", as.double(data), as.integer(n.site),
                       as.integer(n.obs), as.integer(block.size), as.integer(n.block),
                       concProb = double(n.pairs), NAOK = TRUE,
                       PACKAGE = "SpatialExtremes")$concProb

    else
        concProb <- .C("empiricalBootConcProb", as.double(data), as.integer(n.site),
                       as.integer(n.obs), as.integer(block.size),
                       concProb = double(n.pairs), NAOK = TRUE,
                       PACKAGE = "SpatialExtremes")$concProb

    concProb <- pmax(0, concProb)

    if (!missing(n.bins)){
        bins <- c(0, quantile(dist, 1:n.bins/(n.bins+1)), max(dist))
        concProbBinned <- rep(NA, length = n.bins + 1)

        for (k in 1:(n.bins + 1)){
            idx <- which((dist <= bins[k+1]) & (dist > bins[k]))

            if (length(idx)>0)
                concProbBinned[k] <- mean(concProb[idx])
        }

        concProb <- concProbBinned
        dist <- (bins[-1] + bins[-(n.bins+2)])/2
    }

    if (plot){
        if (add)
            points(dist, concProb, col = col[1], ...)

        else
            plot(dist, concProb, ylim = ylim, ylab = expression(p(h)),##expression(p(x[1],x[2])),
                 xlab = expression(h),##xlab = expression(group("||", x[1]-x[2], "||")),
                 xlim = xlim, col = col[1], ...)
        
        if (!missing(fitted)){
            ## Plot the theoretical extremal concurence probability function
            curve(conc.prob.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
        }
    }
        
    return(invisible(list(dist = dist, concProb = concProb)))
}
