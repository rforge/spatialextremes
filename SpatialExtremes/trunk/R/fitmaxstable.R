fitmaxstab <- function(data, coord, cov.mod = c("gauss", "whitmat", "cauchy", "powexp"),
                       loc.form, scale.form, shape.form, fit.marge = TRUE,
                       ..., warn.inf = TRUE, method = "BFGS",
                       std.err.type = "none", corr = FALSE){

  if (missing(loc.form) && missing(scale.form) && missing(shape.form))
    reg.mod <- "full"


  if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form))
    reg.mod <- "spatgev"

  flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)

  if (!(flag %in% c(0, 3)))
    stop("if one formula is given for the GEV parameters, then it should
be given for *ALL* GEV parameters")

  if (cov.mod == "gauss")
    fitted <- switch(reg.mod, "full" = smithfull(data, coord, ..., fit.marge = fit.marge,
                                warn.inf = warn.inf, method = method, std.err.type =
                                std.err.type, corr = corr),
                     "spatgev" = smithform(data, coord, ..., loc.form = loc.form, scale.form = scale.form,
                       shape.form = shape.form, fit.marge = fit.marge,
                       warn.inf = warn.inf, method = method, std.err.type =
                       std.err.type, corr = corr))


  else
    fitted <- switch(reg.mod, "full" = schlatherfull(data, coord, cov.mod = cov.mod,
                                ..., fit.marge = fit.marge, warn.inf = warn.inf,
                                method = method, std.err.type = std.err.type, corr = corr),
                     "spatgev" = schlatherform(data, coord, cov.mod = cov.mod, ...,
                       loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                       fit.marge = fit.marge, warn.inf = warn.inf, method = method, std.err.type =
                       std.err.type, corr = corr))

  return(fitted)
}
