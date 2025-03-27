
#-----------------------------------------------------------------------------------#

#------------------   Gradient Forest Bioregionalisation ---------------------------#

# Fabrice Stephenson & Matt Bennion
# April-May 2024

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
# 2. Hierarchical classification at different levels
# 3. Export classification as rasters
# 4. Analysis of species data

#-----------------------------------------------------------------------------------#
## 1. Load turnover and cluster
#-----------------------------------------------------------------------------------#


setwd("PATH TO FOLDER")
load("Mean_pred.Rdata")
# First -- CLUSTERING 

# Choose number of clusters
n_k <- 50

# CLASSIFY 
StartTime <- Sys.time()
ClaraClassification <- clara(x = Fred.mean,    
                                   k = n_k,
                                   samples = 50,
                                   metric = 'manhattan',   
                                   # Manhattan distance appropriate because transformation 
                                   # -- takes care of scaling differences
                                   trace = 2,
                                   pamLike = T)

EndTime <- Sys.time()
print(EndTime - StartTime)


#-----------------------------------------------------------------------------------#
## 2. Hierarchical classification at different levels
#-----------------------------------------------------------------------------------#

# Now create a 50 row dataset summarising the transformed envrionmental attributes for the groups from clara
MedoidMeans <- matrix(0, nrow = n_k, ncol = length(imp.vars))

dimnames(MedoidMeans)[[1]] <- paste('Grp_',c(1:n_k),sep='')
dimnames(MedoidMeans)[[2]] <- imp.vars

for (i in c(1:length(imp.vars))) MedoidMeans[,i] <- tapply(Fred.mean[,imp.vars[[i]]],
                                                              ClaraClassification[[4]],mean)

# and apply a hierarchical classification to it using agnes
MedoidAgnesClassification <- agnes(MedoidMeans, 
                                      metric = 'manhattan',
                                      method = 'gaverage',
                                      par.method = -0.1)


# now reduce to a smaller number of groups based on the clara results
ClaraGroupExpansion <- cutree(MedoidAgnesClassification,1:n_k)

# head(ClaraGroupExpansion)
i <- match(ClaraClassification$clustering,ClaraGroupExpansion[,n_k])
# summary(i)


Group.Num <- seq(3,n_k,1)
Group.Name <- paste("Grp_",Group.Num, sep = "")

for (j in 1:length(Group.Num)){
  index <- j + 10
  ClaraClassification[[index]] <- ClaraGroupExpansion[i,Group.Num[j]]
  names(ClaraClassification)[index] <- Group.Name[j]
}



#-----------------------------------------------------------------------------------#
## 3. Export classification as rasters
#-----------------------------------------------------------------------------------#


Group.Num <- seq(3,n_k,1)
Group.Name <- paste("Grp_",Group.Num, sep = "")


for (k in 1:length(Group.Num)){
  Means <- matrix(0, nrow = Group.Num[k], ncol = length(imp.vars))
  dimnames(Means)[[1]] <- paste('Grp_', c(1:Group.Num[k]),sep='')
  dimnames(Means)[[2]] <- imp.vars
  
  for (i in c(1:length(imp.vars))) Means[,i] <- tapply(Fred.mean[,imp.vars[[i]]],
                                                       ClaraClassification[[Group.Name[k]]],mean)
  ClusterPCA <- prcomp(Means)
  
  # set up colours using the same PCA space
  a1 <- ClusterPCA$x[, 1]
  a2 <- ClusterPCA$x[, 2]
  a3 <- ClusterPCA$x[, 3]
  r <- a1 + a2
  g <- -a2
  b <- a3 + a2 - a1
  r <- (r - min(r))/(max(r) - min(r)) * 255
  g <- (g - min(g))/(max(g) - min(g)) * 255
  b <- (b - min(b))/(max(b) - min(b)) * 255
  
  ##### Plot PCA ####
  nvs <- dim(ClusterPCA$rotation)[1]
  vec <- imp.vars#c("Mud", "POC", "BPI_broad", "Dissox", "Gravel")  ### How did you select these? Trial and error?
  lv <- length(vec)
  vind <- rownames(ClusterPCA$rotation) %in% vec
  
  scal <- 15
  xrng <- range(ClusterPCA$x[, 1], ClusterPCA$rotation[, 1]/scal) * + 1.1
  yrng <- range(ClusterPCA$x[, 2], ClusterPCA$rotation[, 2]/scal) * + 1.1
  
  variableCEX <- (as.numeric(table(ClaraClassification[[Group.Name[k]]]))^0.110) * 2
  
  # save as EMF
  dir.create(paste(home, "/", Group.Name[k],sep =""))
  setwd(paste(home, "/", Group.Name[k],sep =""))
  jpeg(filename =  paste(Group.Name[k], "_EEZ_PCA.jpeg"), width = 15, height = 15, units = 'cm', quality = 100, bg = "white", res = 500)
  plot((ClusterPCA$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = variableCEX, 
       col = rgb(r, g, b, max = 255), asp = 1)
  arrows(rep(0, lv), rep(0, lv), ClusterPCA$rotation[vec,1]/scal, ClusterPCA$rotation[vec, 2]/scal, length = 0.06210)
  jit <- 0.00110
  text(ClusterPCA$rotation[vec, 1]/scal + jit * sign(ClusterPCA$rotation[vec, 1]), 
       ClusterPCA$rotation[vec, 2]/scal + jit * sign(ClusterPCA$rotation[vec, 2]), labels = vec)
  text(ClusterPCA$x[,1], ClusterPCA$x[,2], seq(1,Group.Num[k],1),
       cex = variableCEX/7, #0.8,
       pos = 4,
       adj = c(0.2,-0.1))
  dev.off() # finishes the plot and saves
  
  ##### Plot map #####
  rFull <- r[ClaraClassification[[Group.Name[k]]]] #clustering]
  gFull <- g[ClaraClassification[[Group.Name[k]]]]
  bFull <- b[ClaraClassification[[Group.Name[k]]]]
  
  jpeg(filename = paste(Group.Name[k], "_Spatial.jpeg", sep = ""), width = 11000, height = 11000, quality = 100, bg = "white", res = 1000)
  par(mar=c(2,2,2,2))
  plot(Pred_1km[,1:2],
       cex = 0.7, 
       col = rgb(rFull, gFull, bFull, max = 255),
       pch = '.',
       asp = 1)
  title(paste(Group.Name[k]))
  dev.off()
  
  #### Create raster ####
  # create a dataframe with spatial information for th classification 
  CMB_DF <- cbind(Pred_1km[,c("x", "y")],ClaraClassification[[Group.Name[k]]])
 
  # export as raster
  CMB_DF.R <- rasterFromXYZ(data.frame(x = CMB_DF[,1],
                                           y = CMB_DF[,2],
                                           z = CMB_DF[,3]),
                                           crs = merc)
  
  writeRaster(CMB_DF.R, filename= paste(Group.Name[k], ".tif", sep = ""), 
              format = "GTiff", 
              overwrite = TRUE)
  
  
  #### Tag speies matrix with class information ####
  class <- as.data.frame(raster::extract(CMB_DF.R, sp.mat[,c("X","Y")]))
  colnames(class) <- Group.Name[k]
  sp.mat <- cbind(sp.mat, class)
}

save(sp.mats, file = "class_tagged.Rdata")


#-----------------------------------------------------------------------------------#
## 4. Analysis of species data
#-----------------------------------------------------------------------------------#

setwd("PATH TO FOLDER")
source("./Pairwise_adonis.R")
require(vegan); require(ecodist); require(parallel)

sp.mat.tag <- sp.mat

Summ.table <- data.frame(Num.Classes =integer(), Prop.total.class=double(),Prop.inter.class.diff=double(),
                         R = double(), p.val = double(), stringsAsFactors=FALSE)


for (j in 1:length(Group.Num)){
  
  grps.sum <- as.data.frame(table(sp.mat.tag[,Group.Name[j]]))
  grps.sum <- grps.sum[grps.sum[,2] > 4,1]
  
  sp.mat.tag.cut <- sp.mat.tag[sp.mat.tag[,Group.Name[j]] %in% grps.sum,]
  
  sp.mat.tag.spe <- sp.mat.tag.cut[,c(21:(ncol(sp.mat.tag.cut)-length(Group.Num)))]   ### Remove Aciniaria? 
  
  if(no.act == T){sp.mat.tag.spe <- sp.mat.tag.spe[,2:17]} # Y/N remove actiniaria - 4385 locations
  
  for (i in 1:ncol(sp.mat.tag.spe)){sp.mat.tag.spe[,i] <- as.numeric(sp.mat.tag.spe[,i])}
  sp.mat.tag.spe[sp.mat.tag.spe==1] <- 0
  sp.mat.tag.spe[sp.mat.tag.spe==2] <- 1
  sp.mat.tag.spe <- sp.mat.tag.spe[,colSums(sp.mat.tag.spe[,1:length(sp.mat.tag.spe)]) > 0]
  sp.mat.tag.spe <- sp.mat.tag.spe[rowSums(sp.mat.tag.spe[1:nrow(sp.mat.tag.spe),]) > 0,]
  
  sp.mat.tag <- sp.mat.tag.cut[,c((ncol(sp.mat.tag.cut)-c(length(Group.Num)-1)):ncol(sp.mat.tag.cut))] ### CHECK THIS
  sp.mat.tag <- sp.mat.tag[rownames(sp.mat.tag) %in% rownames(sp.mat.tag.spe), ]
  
  sp.mat.tag.dist <- vegdist(sp.mat.tag.spe, distance="jaccard", binary = T)

  
  # ANOSIM
  perm <- anosim(sp.mat.tag.dist, sp.mat.tag[,Group.Name[j]], permutations = 100) # added distance
  perm$statistic
  
  Summ.table[j,1] <- Group.Num[j]
  Summ.table[j,2] <- length(grps.sum)/Group.Num[j]
  Summ.table[j,4] <- perm$statistic # Global R
  Summ.table[j,5] <- perm$signif

  print(paste("Iteration finished for ", Group.Num[j], sep =""))
}


write.csv(Summ.table, file = "ANOSIM_R_class_scores.csv")


## Use ggplot2 to visualize the Global R from ANOSIM
library(ggplot2)
df <- data.frame(Summ.table[,1], Summ.table[,4])
p1 <- ggplot(data=df, aes(x=df$Summ.table...1., y=df$Summ.table...4.)) +
  geom_line(color = "#1c9099", linewidth = 1)+
  geom_point(color = "black", size = 2) +
  #geom_vline(xintercept = 6, linetype = 2) + 
  # geom_vline(xintercept = 12, linetype = 2) + 
  # geom_vline(xintercept = 20, linetype = 2) + 
  scale_x_continuous(breaks = seq(5, 50, by = 5)) +
  ylim(0.0,0.20) +
  xlab("Classification level") +
  ylab("Global ANOSIM R-value") +
  theme_bw()
p1

ggsave("Global_R_ANOSIM.png", dpi = 300)

## END ##----------------------------------------------------------------------##

