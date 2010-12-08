library(SpatialExtremes)
library(fields)

##The lambda madogram
n.site <- 50
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")
##Simulate a max-stable process - with unit Frechet margins
data <- rmaxstab(40, locations, cov.mod = "whitmat", sill = 1, range = 1,
smooth = 2)
##Compute the lambda-madogram
par(mar = c(0, 0, 0, 0), pty = "s")
lmadogram(data, locations, n.bins = 80, theta = 45, box = FALSE)
dev.copy2pdf(file = "lmadogram.pdf")
system("convert lmadogram.pdf ../images/lmadogram.png; rm lmadogram.pdf")
dev.off()

##The F-madogram
par(mar = c(4, 5, 1, 1))
fmadogram(data, locations)
dev.copy2pdf(file = "fmadogram.pdf")
system("convert fmadogram.pdf ../images/fmadogram.png; rm fmadogram.pdf")
dev.off()

##A gaussian random field
x <- y <- seq(-10, 10, length = 250)
data <- rgp(1, cbind(x, y), "whitmat", sill = 1, range = 3, smooth = 1.5,
            grid = TRUE)
par(mar = c(4, 4, 1, 1))
image(x, y, data, col = tim.colors(64))
dev.copy2pdf(file = "gp.pdf")
system("convert gp.pdf ../images/gp.png; rm gp.pdf")
dev.off()

##A Smith random field
x <- y <- seq(-10, 10, length = 250)
data <- rmaxstab(1, cbind(x, y), "gauss", cov11 = 3, cov12 = 0.75, cov22 = 3.5,
                 grid = TRUE)
par(mar = c(4, 4, 1, 1))
image(x, y, log(data), col = tim.colors(64))
dev.copy2pdf(file = "smith.pdf")
system("convert smith.pdf ../images/smith.png; rm smith.pdf")
dev.off()

##A Schlather random field
x <- y <- seq(-10, 10, length = 250)
data <- rmaxstab(1, cbind(x, y), "whitmat", sill = 1, range = 3,
                 smooth = 0.75, grid = TRUE)
par(mar = c(4, 4, 1, 1))
image(x, y, log(data), col = tim.colors(64))
dev.copy2pdf(file = "schlather.pdf")
system("convert schlather.pdf ../images/schlather.png; rm schlather.pdf")
dev.off()

##One dimensional Kriging
n.site <- 50
n.pred <- 512
x.obs <- runif(n.site, -100, 100)
x.pred <- seq(-100, 100, length = n.pred)
data <- rgp(1, x.obs, "whitmat", sill = 1, range = 10, smooth = 0.75)
krig <- kriging(data, x.obs, x.pred, "whitmat", sill = 1, range = 10,
smooth = 0.75)

par(mar = c(4, 5, 1, 1))
plot(krig$coord, krig$krig.est, type = "l", xlab = "x", ylab =
 expression(hat(Y)(x)))
points(x.obs, data, col = 2, pch = 21, bg = 2)
dev.copy2pdf(file = "kriging1d.pdf")
system("convert kriging1d.pdf ../images/kriging1d.png; rm kriging1d.pdf")
dev.off()

## Two dimensional kriging on a grid
x.obs <- matrix(runif(2 * n.site, -100, 100), ncol = 2)
x <- y <- seq(-100, 100, length = 100)
x.pred <- cbind(x, y)
data <- rgp(1, x.obs, "whitmat", sill = 1, range = 10, smooth = 0.75)
krig <- kriging(data, x.obs, x.pred, "whitmat", sill = 1, range = 10,
                smooth = 0.75, grid = TRUE)

z.lim <- range(c(data, krig$krig.est))
breaks <- seq(z.lim[1], z.lim[2], length = 65)
col <- tim.colors(64)
idx <- as.numeric(cut(data, breaks))

par(mar = c(4, 4, 1, 1))
image(x, y, krig$krig.est, col = col, breaks = breaks)
points(x.obs, bg = col[idx], pch = 21)
dev.copy2pdf(file = "kriging2d.pdf")
system("convert kriging2d.pdf ../images/kriging2d.png; rm kriging2d.pdf")
dev.off()


##Model checking for a fitted max-stable model
n.site <- 20
n.obs <- 50
coord <- matrix(runif(2 * n.site, 0, 10), ncol = 2)
colnames(coord) <- c("lon", "lat")
data <- rmaxstab(n.obs, coord, "powexp", sill = 1, range = 3, smooth =
                 1)

fitted <- fitmaxstab(log(data), coord, "powexp", y ~ 1, y ~ 1, y ~ 1,
                     sill = 1)
plot(fitted)
dev.copy2pdf(file = "checkMaxStab.pdf")
system("convert checkMaxStab.pdf ../images/checkMaxStab.png; rm checkMaxStab.pdf")
dev.off()

## One dimensional conditional simulation
n.site <- 50
n.sim <- 512
x.obs <- runif(n.site, -100, 100)
x.sim <- seq(-100, 100, length = n.sim)

data <- rgp(1, x.obs, "whitmat", sill = 1, range = 8, smooth = 1.75)
sim <- condrgp(5, x.sim, x.obs, data, "whitmat", sill = 1, range =
               8, smooth = 1.5)

par(mar = c(4, 5, 1, 1))
matplot(x.sim, t(sim$cond.sim),  type = "l", lty = 1, xlab = "x", ylab =
expression(Y[cond](x)))
points(x.obs, data, pch = 21, bg = 1)

dev.copy2pdf(file = "condsim1d.pdf")
system("convert condsim1d.pdf ../images/condsim1d.png; rm condsim1d.pdf")
dev.off()

##Two dimensional simulations
x.obs <- matrix(runif(2 * n.site, -100, 100), ncol = 2)
x <- y <- seq(-100, 100, length = 250)
x.sim <- cbind(x, y)
data <- rgp(1, x.obs, "whitmat", sill = 1, range = 50, smooth = 0.75)

sim <- condrgp(1, x.sim, x.obs, data, "whitmat", sill = 1, range = 50,
               smooth = 0.75, grid = TRUE, control = list(nlines = 2000))

z.lim <- range(c(sim$cond.sim, data, krig$krig.est))
breaks <- seq(z.lim[1], z.lim[2], length = 65)

col <- tim.colors(64)
idx <- as.numeric(cut(data, breaks))

par(mar = c(4, 4, 1, 1))
image(x, y, sim$cond.sim, col = col, breaks = breaks)
points(x.obs, bg = col[idx], pch = 21)


dev.copy2pdf(file = "condsim2d.pdf")
system("convert condsim2d.pdf ../images/condsim2d.png; rm condsim2d.pdf")
dev.off()
