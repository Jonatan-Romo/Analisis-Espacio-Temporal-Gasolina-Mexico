
remove(list = ls())

library(maptools)
library(spdep)
library(leaflet)
library(RColorBrewer)
library(plotly)
library(gstat)
library(rgdal)
library(geoR)
library(ggplot2)
library(spatstat)
library(gapminder)
library(dplyr)
library(plyr)
library(INLA)
library(spacetime)
library(FRK)
library(sp)
library(SpatialEpi)
library(rgeos)
library(MASS)
library(EnvStats)
library(zoo)
library(RMTstat)

setwd("~/Maestr?a en Computo Estad?stico/Tercer semestre/Econometr?a/Proyecto/basesdedatos/")
#Generico <- read.csv("Gen?rico Renta de Vivienda.csv")
Generico <- read.csv("Genérico Gasolina de Alto Octanaje.csv")
times <- nrow(Generico)/46
Generico <- Generico[,c(1,2,3,4,8)]

index <- rep(1:46,times)
std <- vector()
i <- 1
for (i in 1:46) {
   std[i] <- var(Generico$Precio.promedio[index==i])
 }
Generico$std <- rep(std,times)
Generico <- within(Generico,
               {time = as.Date(as.yearmon(paste(Generico$A?o,Generico$Mes,sep="-")))}) # create Date field
proy <- SpatialPoints(Generico[,c(3,4)], CRS("+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
coords <- spTransform(proy,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
Generico[,c(3,4)] <- coords@coords

ruta <- "~/Maestr?a en Computo Estad?stico/Tercer semestre/Econometr?a/Proyecto/SHP"
Local <- readOGR(dsn=ruta, layer = ogrListLayers(ruta)[1])

ruta <- "~/Maestr?a en Computo Estad?stico/Tercer semestre/Econometr?a/Proyecto/MexicoEstados"
Mexico.borders <- readOGR(dsn=ruta, layer = ogrListLayers(ruta)[1])
Mexico.borders.Tr <- spTransform(Mexico.borders,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
Union <- spTransform(unionSpatialPolygons(Mexico.borders,rep(1,32)),"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

plot(Mexico.borders.Tr, axes=TRUE)
plot(coords, add = TRUE, pch = 16, col = "blue")

STObj <- stConstruct(x = Generico, # dataset
                     space = c("lon","lat"), # spatial fields
                     time="time", # time field
                     interval=TRUE) # time reflects an interval

grid_BAUs <- auto_BAUs(manifold=STplane(), # spatio-temporal process on the plane
                       data=STObj, # data
                       cellsize = c(0.9,0.9,1), # BAU cell size
                       type="grid", # grid or hex?
                       convex=-0.15, # border buffer factor
                       nonconvex_hull=TRUE, # convex hull
                       tunit="months") # time unit
grid_BAUs$fs = 1 # fine-scale variation
grid_BAUsR <- grid_BAUs

plot(grid_BAUs@sp@coords, pch = 16, col = "blue")
plot(Mexico.borders.Tr, add=TRUE)

G_spatial <- auto_basis(manifold = plane(), # spatial functions on the plane
                        data=as(STObj,"Spatial"), # remove the temporal dimension
                        nres = 2, # 2 resolutions
                        type = "Gaussian", # bisquare basis functions
                        regular = 0,
                        max_basis = 1000,
                        verbose = TRUE) # regular basis functions

G_temporal <- local_basis(manifold = real_line(), # functions on the real line
                          loc=matrix(1:10,10,1),
                          scale=rep(2,10),
                          type="Gaussian") # scales of functions

basis_s_plot <- show_basis(G_spatial) + xlab("lon") + ylab("lat")
basis_t_plot <- show_basis(G_temporal) + xlab("time index") + ylab(expression(phi(t)))

G_spacetime <- TensorP(G_spatial,G_temporal) # take the tensor product

f <-Precio.promedio ~ 1 + lat + lon # fixed effects part
Frk <- FRK(f = f, # formula
         data = list(STObj), # spatio-temporal object
         basis = G_spacetime, # space-time basis functions
         BAUs = grid_BAUs, # space-time BAUs
         est_error = FALSE, #  measurement error                     
         average_in_BAU = TRUE, # average data that fall inside BAUs
         tol = 0.01, n_EM = 1)

Pred.Grid <- predict(Frk, newdata = Prueba) 
Pred.Grid
stplot(Pred.Grid)
Pred.Grid.Pr <- Pred.Grid 
Pred.Grid.Pr@data <- Pred.Grid.Pr@data[,-c(1,3,4)]
Pred.Grid.Pr@data <- Pred.Grid.Pr@data[,c(2,3,1)]

Mex <- over(SpatialPoints(Pred.Grid.Pr@sp@coords, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")),Union)
Pred.Grid.Pr@sp@coords <- Pred.Grid.Pr@sp@coords[!is.na(Mex),]
Pred.Grid.Pr@sp@data <- Pred.Grid.Pr@sp@data[!is.na(Mex),]
Pred.Grid.Pr@sp@grid.index <- Pred.Grid.Pr@sp@grid.index[!is.na(Mex)]
Pred.Grid.Pr@data <- Pred.Grid.Pr@data[rep(!is.na(Mex),(times+1)),]

stplot(Pred.Grid.Pr)

st.pr <- stplot(Pred.Grid.Pr, main="Precios en M?xico",
       scales=list(draw=F),
       sp.layout = list(list("sp.polygons", Union, first=FALSE,col="black"),
       list("sp.points", coords, col=gray(0.25), pch=3, cex=.5)))
plot(st.pr)

Pred.Grid.Pr.Se <- Pred.Grid.Pr
Pred.Grid.Pr.Se@data <- Pred.Grid.Pr.Se@data[,c(2,1)]

st.pr.se <- stplot(Pred.Grid.Pr.Se,
       scales=list(draw=F),
       sp.layout = list(list("sp.polygons", Union, first=FALSE,col="black"),
                        list("sp.points", coords, col=gray(0.25), pch=3, cex=.5)))
plot(st.pr.se)

Pred <- as(Pred.Grid.Pr,"data.frame")
u.coords <- unique(Generico[,c(3,4)])
u.grid <- unique(Pred[,c(1,2)])
Dist.Matrix <- distR(u.coords,u.grid)
coords.rel <- sapply(s <- 1:nrow(u.coords), function(s) { which.min(Dist.Matrix[s,]) })
near.grid <- u.grid[coords.rel,]

Predicho <- vector()
for(s in 2:20){
  Predicho <- c(Predicho,Pred[Pred$t == s,9][coords.rel])
}
se <- (Generico[,5] - Predicho)^2
se.matrix <- matrix(se, nrow = 19, byrow = TRUE)
MSE <- data.frame(x= coords@coords[1:46,1], y = coords@coords[1:46,2], mse = colMeans(se.matrix))
mean(se)

# add to data a new column termed "id" composed of the rownames of data
Mexico.borders.Tr@data$id <- rownames(Mexico.borders.Tr@data)

# create a data.frame from our spatial object
MexicoPoints <- fortify(Mexico.borders.Tr, region = "id")

# merge the "fortified" data with the data from our spatial object
MexicoDF <- merge(MexicoPoints, Mexico.borders.Tr@data, by = "id")
head(MexicoDF)

ggMexicoDF <- ggplot() +
  geom_polygon(data = MexicoDF, aes(x=long, y=lat, group = group), color = "gray", fill = "white")  +
  #geom_path(color = "white") +
  #scale_fill_hue(l = 40) +
  coord_equal() +
  #theme(legend.position = "none", title = element_blank(),
  #      axis.text = element_blank()) +
  geom_point(data=MSE, # Plot data
             aes(x,y, fill = mse), # Colour <-> log(zinc)
             colour="blue", # point outer colour
             pch=21, size=3) # size of point
print(ggMexicoDF)

error <- (Generico[,5] - Predicho)
error.matrix <- matrix(error, nrow = 19, byrow = TRUE)
  

# Ahora, de la tabla del art?culo original de Stephen W. Looney y Thomas R. Gulledge. Para $n = 46$ y $\alpha = 0.05$ el punto cr?tico de la correlaci?n es $0.975$.
#Siendo $r_Q=0.9861917 > 0.975$, por lo que no rechazamos la hip?tesis de normalidad.
# Los ?ltimos 4 puntos pertenecen a puntos de la frontera... Procedemos a un an?lisis con matrices aleatorias
# Escalamos para el analisis con matrices aleatorias y sus valores propios
Scaled.Error <- scale(error.matrix)
N <- dim(Scaled.Error)[1]
P <- dim(Scaled.Error)[2]

A <- t(Scaled.Error) %*% Scaled.Error / P
# Es b?sicamente una estimaci?n de la varianza,

# Mostramos la distribuci?n de Tracy-Widom con el limite superior de 0.05 de significancia
# para considerar p-value < 0.05
s <- seq(-5,5,0.1)
lim <- qtw(0.05, lower.tail = FALSE)


# armamos el estad?stico y determinamos aquellos valores mayores que el limite

lambda <- eigen(A)$values
mu <- (sqrt(N-0.5)+sqrt(P-0.5))^2
sn <- mu*(1/sqrt(N-0.5)+1/sqrt(P-0.5))^(1/3)
Componentes <- sum((N*lambda-mu)/sn > lim)

# Mostramos el histograma de los valores propios (omitiendo el mayor) con el limite marcado
plot(hist(((N*lambda-mu)/sn),breaks = seq(-3,5,0.5)), xlim = c(-3,5))
abline(v = lim, col = "red")

# Mostramos en rojo cuando la acumulada de TW llega 0.95 y en azul el estad?stico obtenido
plot(s, dtw(s), type = "lines", ylab = "Probabilidad")
abline(v = lim, col = "red")
abline(v = ((N*lambda-mu)/sn)[1], col = "blue")
ptw(((N*lambda-mu)/sn)[1], lower.tail = FALSE)

# Mostramos una proyecci?n de los residuales resaltando en rojo los datos m?s ajenos a la muestra
comp <- prcomp(A)
plot(comp$x[,1:2], pch = 18)
points(comp$x[c(5,6,24,26,40),1:2],col = "red", pch = 18)


# Mostramos en una imagen geogr?fica las localidad mostradas en rojo en la anterior grafica
ggMexicoDFObs <- ggplot() +
  geom_polygon(data = MexicoDF, aes(x=long, y=lat, group = group), color = "gray", fill = "white")  +
  #geom_path(color = "white") +
  #scale_fill_hue(l = 40) +
  coord_equal() +
  #theme(legend.position = "none", title = element_blank(),
  #      axis.text = element_blank()) +
  geom_point(data=MSE[c(5,6,24,26,40),], # Plot data
             aes(x,y, fill = mse), # Colour <-> log(zinc)
             colour="blue", # point outer colour
             pch=21, size=3) # size of point
print(ggMexicoDFObs)



# Observamos que son los datos de la frontera los que son an?malos, y esto es de esperarse
# ya que en la frontera entro m?s competencia.
# Por la liberaci?n y la cercan?a con EEUU. Dado que la gasolina se compra en gran parte a EEUU, quiz? eso influye
