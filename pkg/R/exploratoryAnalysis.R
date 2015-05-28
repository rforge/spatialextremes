symbolplot <- function(data, coord, which = "gev", plot.border = NULL, col = c("#FF000080", "#0000FF80")){

    if (!(which %in% c("gev", "mean", "median")))
        stop("'which' must be one of 'gev', 'mean' or 'median'")


    add <- FALSE
    if (which == "gev"){
        values <- apply(data, 2, gevmle)
        ref.value <- rowMeans(values)
        
        par(mfrow = c(1, 3))

        for (i in 1:3){
            if (!is.null(plot.border)){
                plot.border()
                add <- TRUE
            }

            sign.col <- col[2 - (values[i,] >= ref.value[i])]
            radius <- abs(values[i,] - ref.value[i])
            norm.factor <-  min(dist(coord)) / max(radius)
            radius <- norm.factor * radius
            symbols(coord, circles = radius, add = add, bg = sign.col, inches = FALSE)
        }
    }

    else {
        if (which == "mean")
            values <- colMeans(data, na.rm = TRUE)
        else
            values <- apply(data, 2, median, na.rm = TRUE)

        ref.value <- mean(values, na.rm = TRUE)

        if (!is.null(plot.border)){
            plot.border()
            add <- TRUE
        }

        sign.col <- col[2 - (values >= ref.value)]
        radius <- abs(values - ref.value)
        norm.factor <- min(dist(coord)) / max(radius)
        radius <- norm.factor * radius

        symbols(coord, circles = radius, add = add, bg = sign.col, inches = FALSE)
    }

    par(new = FALSE)
}


