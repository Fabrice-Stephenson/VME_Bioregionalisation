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

# 1. Load packages and biological and environmental data
# 2. Create environmental ranges
# 3. Run bootstrapped gradient forest model
# 4. Post-bootstrap processing
# 5. PCA and generation of raster outputs

# Clustering (script no. 2)

#-----------------------------------------------------------------------------------#
## 1. Load packages
#-----------------------------------------------------------------------------------#

require(extendedForest)
require(gradientForest)
require(cluster)
require(raster)
library(tidyr)
library(dplyr)
library(plyr)
library(foreach)
library(data.table)


## Clear whole environment
##  rm(list = ls()) 


#-----------------------------------------------------------------------------------#
# 1. Load biological and environmental data
#-----------------------------------------------------------------------------------#

# Worms accepted taxa list - created AFTER FID and matrix creation
setwd("PATH TO FOLDER")
taxa_list <- read.csv("Taxa_to_model.csv")
taxa <- taxa_list$taxa

load("env_pred.Rdata") # environmental data for prediction ## Name: Pred_1km
load("species_matrix.Rdata") # biological data matrix ## Name: sp.mat

# Set a CRS
merc <- CRS("+proj=merc +lat_ts=-41 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Create an object with names of the environmental variables to use
imp.vars <- c("BPI_broad", "Slope_SD", "Ruggedness", "Dissox","Aragonite",
              "Temp", "POC", "Depth", "Silicate", "Salinity") 



#-----------------------------------------------------------------------------------#
# 2. Create environmental ranges
#-----------------------------------------------------------------------------------#

# create environmental gradient files for each taxa
# first subset Pred_1km (environmental prediction dataset) by imp.vars

# skip first two columns - XY coords
PredMins <- apply(Pred_1km[,c(3:ncol(Pred_1km))],2,min)
PredMaxs <- apply(Pred_1km[,c(3:ncol(Pred_1km))],2,max)

EnvRanges <- as.data.frame(array(0,c(200,length(imp.vars))))
names(EnvRanges) <- imp.vars

for (i in c(1:length(imp.vars))) {
  EnvRanges[,i] <- seq(PredMins[imp.vars[i]],PredMaxs[imp.vars[i]],length = 200)
}


#-----------------------------------------------------------------------------------#
# 3. Bootstrapped gradient forest model
#-----------------------------------------------------------------------------------#

# Set a path to a folder to save GF model objects
setwd("PATH TO FOLDER")

# Number of bootstraps - 100
n.boot <- 100

# Create saving objects
multiResultClass <- function(ImpVar.fam=NULL,SpeR2=NULL, 
                                 EnvTran=NULL, TurnOv=NULL,
                                 model = NULL ) {me <- list(
                                 ImpVar.fam = ImpVar.fam, SpeR2 = SpeR2, 
                                 EnvTran = EnvTran, TurnOv = TurnOv, model = model)
                                 ## Set the name for the class
                                 class(me) <- append(class(me),"multiResultClass")
                                 return(me)
    }
    
    library(bigstatsr)
    cl <- parallel::makeCluster(11) # number of cores used
    # showConnections()
    doParallel::registerDoParallel(cl)
    packages <- c("extendedForest", "gradientForest") 
    # need to call packages within each parallel loop
    
    start <- Sys.time() # recording the time
    
GF.boots <- foreach(i = 1:n.boot) %dopar% { 
      
# Gradient forest model
nSites <- nrow(sp.mat)
lev <- floor(log2(nSites * 0.368/2))
# create training dataset

train_ind <- sample(seq_len(nrow(sp.mat)), size = nSites, replace = T) 
DF_train <- sp.mat[train_ind,]

# Model  colnames(DF_train)
lapply(packages, require, character.only = TRUE)
GF_obj <- gradientForest(DF_train,
                                predictor.vars = imp.vars,
                                response.vars = taxa,
                                ntree = 250, 
                                transform = NULL,
                                compact = T,
                                nbin = 401,
                                maxLevel = lev,
                                corr.threshold = 0.5,
                                trace = T)


result <- multiResultClass()

result$model <- GF_obj
# Imp of envir variables (sorted into alphabetical order)
result$ImpVar.fam <- sort(round(GF_obj$overall.imp, 5))
# Species model fits 
result$SpeR2 <- GF_obj$result
# Env transforms
result$EnvTran <- predict(GF_obj, EnvRanges, extrap=F) 
result$TurnOv <- predict(GF_obj, Pred_1km[,imp.vars], extrap=F)

# Return results files
return(result)
}
    
parallel::stopCluster(cl)
    
end <- Sys.time()
end - start
    
setwd("./bootstrap_results")
save(GF.boots, file="GF.boots.Rdata")
#load("GF.boots")

#-----------------------------------------------------------------------------------#
# 4. Post-bootstrap processing
#-----------------------------------------------------------------------------------#

## Next section processes bootstraps, produces mean turnover and SD and extracts fits
#list of files
ls <- list("1" = "GF.boots")

for (a in 1:length(ls)){
  # load output from first set of bootstraps
  setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results", sep = ""))
  load(paste(ls[[a]],".Rdata", sep = "")) # GF.boots <- 
  
  # Important variables
  boot_array_ImpVar <- array(0, dim = c(length(imp.vars), n.boot))
  rownames(boot_array_ImpVar) <- imp.vars
  
  for (i in 1:n.boot){
  boot_array_ImpVar[,i] <- GF.boots[[i]]$ImpVar.fam[imp.vars]
  }
  # loop through all versions to save full set of boot array
  if (a == 1){boot_array_ImpVar_F <- boot_array_ImpVar
  } else {boot_array_ImpVar_F <- cbind(boot_array_ImpVar_F, boot_array_ImpVar)}
  
  # save final object 
  if (a == length(ls)){preds_influences <-t(apply(boot_array_ImpVar, 1, function(x) c(Mean = mean(x), SD = sd(x))))   #Calculate mean and standard error of the relative influecne of each preds)
  setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
  write.csv(preds_influences, file="preds.influences_genus_level.csv")}
  
  # SPECIES ## JESUS maybe fix the colnames bit
  boot_array_SpeR2 <- array(0, dim = c(length(taxa), n.boot)) # refer back to the  full dataset
  rownames(boot_array_SpeR2) <- sort(taxa)
  
  for (i in 1:n.boot){
  boot_array_SpeR2[,i] <- GF.boots[[i]]$SpeR2[taxa]
  }
  
  # Number of species modelled, min, max and mean fit per GF boot
  SpeR2_sum <- array(0, dim = c(n.boot, 4)) # refer back to the  full dataset
  colnames(SpeR2_sum) <- c("length", "min", "mean","max")
  
  for (i in 1:n.boot){
  SpeR2_sum[i,1] <- length(boot_array_SpeR2[,i]) - sum(is.na(boot_array_SpeR2[,i]))
  SpeR2_sum[i,2] <- min(boot_array_SpeR2[,i], na.rm = T)
  SpeR2_sum[i,3] <- mean(boot_array_SpeR2[,i], na.rm = T)
  SpeR2_sum[i,4] <- max(boot_array_SpeR2[,i], na.rm = T)
  }
  
  # loop through all versions to save full set of boot array
  if (a == 1){SpeR2_sum.fam_F <- SpeR2_sum
  } else {SpeR2_sum.fam_F <- rbind(SpeR2_sum.fam_F, SpeR2_sum)}
  # save final object 
  if (a == length(ls)){SpeR2_sum.fam_F <- round(apply(SpeR2_sum.fam_F, 2, function(x) c(Mean = mean(x), SD = sd(x))),2)
  setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
  write.csv(SpeR2_sum.fam_F, file="Sp_meanR2_genus_level.csv")}
  
  # species model fit and variability
  # loop through all versions to save full set of boot array
  boot_array_SpeR2[is.na(boot_array_SpeR2)] <- 0
  if (a == 1){boot_array_SpeR2_F <- boot_array_SpeR2
  } else {boot_array_SpeR2_F <- cbind(boot_array_SpeR2_F, boot_array_SpeR2)}
  
  # save final object 
  if (a == length(ls)){speciesR2 <- t(apply(boot_array_SpeR2_F, 1, function(x) c(Mean.R2 = mean(x), SD = sd(x))))
  setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
  write.csv(speciesR2, file="Sp_R2.csv")}
  
  # PREDICTED ENVIRONMENTAL TRANSFORMS
  boot_array_EnvTran <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot)) # EnvRanges NEEDED
  dimnames(boot_array_EnvTran)[[2]] <- imp.vars
  dimnames(boot_array_EnvTran)[[3]] <- paste('Rep_',seq(1,n.boot),sep="")
  
  for (i in 1:n.boot){
    
  for (j in 1:length(EnvRanges)){
  boot_array_EnvTran[,j,i] <- GF.boots[[i]]$EnvTran[,j]
  }
  }
  
  if (a == 1){boot_array_EnvTran <- boot_array_EnvTran
  } else {boot_array_EnvTran <- abind::abind(boot_array_EnvTran, boot_array_EnvTran, along = 3)
  }
  # save final object 
  if (a == length(ls)){setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
  save(boot_array_EnvTran, file="boot_array_EnvTran.Rdata")}
  
  # FRED
  Fred <- array(0,c(length(Pred_1km[[1]]),length(EnvRanges))) # EnvRanges
  dimnames(Fred)[[2]] <- imp.vars
  
  # CALCULATE MEAN 
  for (i in 1:n.boot){Fred <- Fred + GF.boots[[i]]$TurnOv} 
  Fred.1 <- Fred/n.boot
  
  if (a == 1){Fred.F <- Fred.1} else {Fred.F <- Fred.F+Fred.1}
  # save final object 
  if (a == length(ls)){Fred.mean <- Fred.F / a
  setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
  save(Fred.mean, file="Mean_pred.Rdata")}
  
  # FRED 2 - calculate SD 
  Fred.2 <- array(0,c(length(Pred_1km[[1]]),length(EnvRanges)))
  dimnames(Fred.2)[[2]] <- imp.vars
  
  # Calculate the mean 
  for (i in 1:n.boot){
    if (i == 1){
      Fred.2 <- (Fred.1 - tmp3[[1]]$TurnOv)^2
    } else {Fred.2 <- Fred.2 + (Fred.1 - tmp3[[i]]$TurnOv)^2
    }
  }
  
  Fred.2.sq <- sqrt(Fred.2)
  
  if (a == 1){Fred.2.f <- Fred.2.sq} else {Fred.2.f <- Fred.2.f + Fred.2.sq}
  # save final object 
  if (a == length(ls)){Fred.2.f_mean <- Fred.2.f / a
  Fred.2.f_mean <- as.vector(apply(Fred.2.f_mean, 1, mean)) 
  setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
  save(Fred.2.f_mean, file="SD_pred.Rdata")} 
  
}


#-----------------------------------------------------------------------------------#
# 5. PCA and generation of raster outputs
#-----------------------------------------------------------------------------------#

setwd("PATH TO FOLDER")

setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/bootstrap_results/after_boots", sep = ""))
load("Mean_pred.Rdata")

# SETUP PCA 
PCs <- prcomp(na.omit(Fred.mean))
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
# Setup plotting parameters
nvs <- dim(PCs$rotation)[1]
vec <- imp.vars
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <- 12 # Bigger number zoom in
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * + 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * + 1.1

# Plot in the PCA space
# windows()
jpeg(filename = "PCA_VME_GF.jpeg", width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot((PCs$x[, c(1,2)]), xlim = xrng, ylim = yrng, pch = ".", cex = 2, col = rgb(r, g, b, max = 255), asp = 1)
#points(PCs$rotation[!vind, c(1,2)]/scal, pch = "+")
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec, 2]/scal, length = 0.0625)
jit <- 0.0015
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec, 1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec, 2]), labels = vec)
title('VME GF model')
dev.off()


# Now plot the PCA in geographic space via RGB
jpeg(filename = "Spatial_Compositional_Turnover.jpeg", 
     width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
par(mar=c(2,2,2,2))
plot(Pred_1km[, c("x", "y")], pch = ".", cex = 1, asp = 1, col = rgb(r, g, b, max = 255))
title('GF Compositional Turnover')
dev.off()

# save files
GF_ras_R <- rasterFromXYZ(data.frame(x = Pred_1km[, c("x")], 
                                     y = Pred_1km[, c("y")], 
                                     z = r),
                          crs = merc)
GF_ras_G <- rasterFromXYZ(data.frame(x = Pred_1km[, c("x")], 
                                     y = Pred_1km[, c("y")], 
                                     z = g),
                          crs = merc)
GF_ras_B <- rasterFromXYZ(data.frame(x = Pred_1km[, c("x")], 
                                     y = Pred_1km[, c("y")],
                                     z = b),
                          crs = merc)

# Write raster tifs - can then be mapped in GIS via RGB
writeRaster(GF_ras_R,"R_mean.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_G,"G_mean.tif","GTiff", overwrite=TRUE)
writeRaster(GF_ras_B,"B_mean.tif","GTiff", overwrite=TRUE)

## End of Script 1
## See Script no. 2 for clustering step
#-----------------------------------------------------------------------------------#
