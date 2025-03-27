#-----------------------------------------------------------------------------------#

#------------------   Gradient Forest Bioregionalisation ---------------------------#

# Fabrice Stephenson & Matt Bennion
# April/May 2024

#-----------------------------------------------------------------------------------#

# This code was produced to develop a VME bioregionalisation for the High
# Seas in the Western part of the South Pacific which is described in detail in
# Bennion et al., (in review) ‘A bioregionalisation for vulnerable marine ecosystems
# in the South Pacific Ocean’

# The analysis adapts methods and code described in:
# 1. Stephenson et al. (2022). Development of a Seafloor Community
#     Classification for the New Zealand Region Using a Gradient Forest Approach
#     Frontiers in Marine Science 8.
# 2. Stephenson et al. (2023). A seafloor bioregionalisation for New Zealand.
#     Ocean & Coastal Management 242, 106688.
# Code for these is available at: https://github.com/Fabrice-Stephenson

#-----------------------------------------------------------------------------------#

# 1. Load turnover and cluster
# 2. Estimate coverage of environmental space

#-----------------------------------------------------------------------------------#
## 1. Load packages and files and prepare dataframe
#-----------------------------------------------------------------------------------#

require(raster)
require(dismo);require(gbm); require(devEMF)

# set CRS
merc <- CRS("+proj=merc +lat_ts=-41 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# load master predictor stack 
setwd("PATH TO FOLDER")
load("Pred_1km.Rdata") # (1 km resolution environmental dataframe)

## Also load species dataframe
load("DF.Rdata")

#-----------------------------------------------------------------------------------#
## 2. Estimate coverage of environmental space
#-----------------------------------------------------------------------------------#

# coverage by absences - first create raster to extract cell ID - FID
r <- rasterFromXYZ(data.frame(X =Pred_1km[,1], 
                              Y =Pred_1km[,2], 
                              Z =1),
                              crs = crs)

Pred_1km_U <-  Pred_1km
Pred_1km_U$FID <-  cellFromXY(r, Pred_1km_U[,c("x","y")])
Pred_1km_U <- Pred_1km[!Pred_1km_U$FID %in% DF$FID,] # no overlap with presences

set.seed(314)
preddat_EC <- Pred_1km_U[sample(seq_len(nrow(Pred_1km_U)), size = nrow(DF)),]
preddat_EC$pa <- 0 # absences in Env Space
env_vars <- PA[,imp.var]
x_y <- DF[,c("X", "Y")]
PA_bind <- cbind(env_vars, x_y)
colnames(preddat_EC)[1] <- "X"
colnames(preddat_EC)[2] <- "Y"
PA_bind$pa <- 1


DF_ES <- rbind(preddat_EC, PA_bind)


DF_BRT <- gbm.step(data=DF_ES, gbm.x = imp.var, # environmental variable columns (imp.var)
                      gbm.y = c("pa"), # presence-absence column
                      family = "bernoulli", tree.complexity = 2,
                      learning.rate = 0.05, bag.fraction = 0.6, n.folds=10, 
                      max.trees = 5000, plot.main = T, tolerance.method = "fixed",
                      tolerance = 0.01, verbose = T)

pred.map <- predict.gbm(DF_BRT, Pred_1km, 
                        n.trees = DF_BRT$gbm.call$best.trees, type = "response")
# export the map
ES.map.mean <- cbind(Pred_1km[, c("x", "y")],pred.map)

# convert to raster
BRT_ES.mean <- rasterFromXYZ(data.frame(x = ES.map.mean[,1], 
                                        y = ES.map.mean[,2], 
                                        z = ES.map.mean[,3]),
                                        crs = merc) 

setwd("PATH TO FOLDER")
writeRaster(BRT_ES.mean,filename = "Env_coverage.tif", overwrite=T) # write raster file

#-----------------------------------------------------------------------------------#

