# Scott Merrill
# Practicing using spatial data
# https://rpubs.com/nabilabd/118172

rm(list=ls(all=T))

install.packages("gstat")
install.packages("geoR")
install.packages("ncf")

require(datasets)
require(sp)
require(gstat) # install.packages("gstat")
require(geoR) # install.packages("geoR")
require(raster)
library(spdep)
require(ncf) # install.packages("ncf")
require(ggplot2)

library(dplyr) # for "glimpse"
library(scales) # for "comma"
library(magrittr)




# https://rpubs.com/nabilabd/118172
demo(meuse, ask = FALSE, echo = FALSE)
head(meuse)
str(meuse) # problem - creating a SpatialPointsDataFrame
coordinates(meuse) # can you redefine coordinates here? i.e,. coordinates(meuse) = ...
coordinates(meuse) = seq(length.out=length(coordinates(meuse)),0,100) # apparently not :(

? variogram


v = variogram(log(zinc)~elev, meuse) # object is a SpatialPointsDataFrame, format model~regressors (independent variables)
plot(v)


v.fit = fit.variogram(v, vgm(1, c("Exp", "Mat", "Sph", "Gau"), 700, 1))
# original model example uses only "Sph" above
# ? vgm
# Need to remember how to fit the best model, Spherical, Gaussian, Exponential or Matern "Mat"
# or use fit.variogram(v, vgm(c("Exp", "Mat", "Sph"))) and let the computer figure it out
v.fit

lzn.vgm <- variogram(log(zinc)~1, meuse) # calculates sample variogram values. Use on Residuals instead of log(zinc)
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1)) # fit model

plot(lzn.vgm, lzn.fit)

# load spatial domain to interpolate over
data("meuse.grid")

# to compare, recall the bubble plot above; those points were what there were values for. this is much more sparse
plot1 <- meuse  %>% as.data.frame %>%
  ggplot(aes(x, y)) + geom_point(size=1) + coord_equal() +
  ggtitle("Points with measurements")

# this is clearly gridded over the region of interest
plot2 <- meuse.grid %>% as.data.frame %>%
  ggplot(aes(x, y)) + geom_point(size=1) + coord_equal() +
  ggtitle("Points at which to estimate")

lzn.kriged <- krige(log(zinc) ~ 1, meuse, meuse.grid, model=v.fit)
plot(lzn.kriged)
spplot(lzn.kriged)


library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)

lzn.kriged2 = lzn.kriged %>% as.data.frame %>%
  ggplot(aes(x=x, y=y)) + geom_tile(aes(fill=var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high="red") +
  scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma) +
  theme_bw()

meuse %>% as.data.frame %>%
  ggplot(aes(x, y)) + geom_point(aes(size=zinc), color="blue", alpha=3/4) +
  ggtitle("Zinc Concentration (ppm)") + coord_equal() + theme_bw()

set = list(gls=1)

g = gstat(NULL, "log-zinc", log(zinc)~x+y, meuse, model=v.fit, set = set)
variogram(g)

if (require(rgdal)) {
  proj4string(meuse) = CRS("+init=epsg:28992")
  meuse.ll = spTransform(meuse, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
  # variogram of unprojected data, using great-circle distances, returning km as units
  variogram(log(zinc) ~ 1, meuse.ll)
}
