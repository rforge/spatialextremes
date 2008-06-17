fitmaxstab <- function(data, coord, cov.mod = c("gauss", "whitmat", "cauchy", "powexp"),
                       loc.form, scale.form, shape.form, fit.marge = TRUE,
                       marg.cov = NULL, ..., warn.inf = TRUE, method = "BFGS",
                       std.err.type = "score", corr = FALSE){

  if (missing(loc.form) && missing(scale.form) && missing(shape.form))
    reg.mod <- "full"


  if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form)){
    reg.mod <- "spatgev"

    if ((class(loc.form) != "formula") || (class(scale.form) != "formula") ||
        (class(shape.form) != "formula"))
      stop("``loc.form'', ``scale.form'' and ``shape.form'' must be valid R formulas")
  }

  flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)

  if (!(flag %in% c(0, 3)))
    stop("if one formula is given for the GEV parameters, then it should
be given for *ALL* GEV parameters")

  if (cov.mod == "gauss")
    fitted <- switch(reg.mod, "full" = smithfull(data, coord, ..., fit.marge = fit.marge,
                                warn.inf = warn.inf, method = method, std.err.type =
                                std.err.type, corr = corr),
                     "spatgev" = smithform(data, coord, ..., loc.form = loc.form, scale.form = scale.form,
                       shape.form = shape.form, fit.marge = fit.marge, marg.cov = marg.cov,
                       warn.inf = warn.inf, method = method, std.err.type =
                       std.err.type, corr = corr))


  else
    fitted <- switch(reg.mod, "full" = schlatherfull(data, coord, cov.mod = cov.mod,
                                ..., fit.marge = fit.marge, warn.inf = warn.inf,
                                method = method, std.err.type = std.err.type, corr = corr),
                     "spatgev" = schlatherform(data, coord, cov.mod = cov.mod, ...,
                       loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                       fit.marge = fit.marge, marg.cov = marg.cov, warn.inf = warn.inf,
                       method = method, std.err.type = std.err.type, corr = corr))

  return(fitted)
}
