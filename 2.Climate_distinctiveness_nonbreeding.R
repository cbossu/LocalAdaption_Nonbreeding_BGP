library(tidyverse)
library(rworldmap)
library(sf)
library(raster)
library(terra)
library(fields)
library(ggpubr)
library(ebirdst)
library(tidyterra)
library(MASS)
library(pracma)
library(ggnewscale)
library(gridExtra)
library(xlsx)

setwd("~/Rproj_git/LocalAdaption_Nonbreeding_BGP/")
## This is the metadata obtained from the previously conducted genoscape analyses. It has 5 columns: species, season, population, longitude, and latitude. Each row is an individual sample.
data_for_analysis <- read.csv("input/distinctiveness/data_for_analysis_new_wintering.3species.csv")
species.names <- unique(data_for_analysis$species)

## Map of the Americas
newmap <- getMap(resolution = "low")
newmap <- spTransform(newmap, '+proj=longlat +datum=WGS84')
newmap@data$world <- rep(1,length(newmap@data$SOVEREIGNT))
newmap <- vect(newmap)
newmap <- terra::buffer(newmap, 0)
newmap <- terra::aggregate(newmap)
newmap <- terra::crop(newmap, ext(-180,-30,-60,90))

## Load ecoregion polygons (obtained from Dinerstein et al. 2017 Biosciences — freely available at https://ecoregions.appspot.com/)) and too large to upload here
ecoregions <- terra::vect("input/distinctiveness/resources/Ecoregions2017/Ecoregions2017.shp")
crs(ecoregions) <- "+proj=longlat +datum=WGS84"
ecoregions <- terra::crop(ecoregions, newmap)
ecoregions <- ecoregions[-96,] # remove an ecoregion whose polygon is causing problems with spatial manipulations


## Download ebird seasonal relative abundance surfaces for all species (freely available)
## Need a key to access eBird data
set_ebirdst_access_key("XXXXX",overwrite=T)
species.names.ebird <- c("wilfly", "swathr", "amered")
sp_path <- "input/distinctiveness/resources/ebird_st/2023"
for(i in 1:length(species.names.ebird)){
  ebirdst::ebirdst_download_status(species.names.ebird[i], path="input/distinctiveness/resources/ebird_st")
}

## Calculate the species relative abundance in each ecoregion and for each species
ecoregions_abund_presence <- list()
for(i in 1:length(species.names)){
  data_for_analysis2 <- data_for_analysis %>% filter(species == species.names[i])
  pops <- unique(data_for_analysis2$population)
  
  # Load eBird seasonal abundance surfaces
  abunds <- terra::rast(paste0(sp_path, "/", species.names.ebird[i], "/seasonal/", species.names.ebird[i], "_abundance_seasonal_mean_9km_2023.tif"))
  abunds <- terra::project(abunds, '+proj=longlat +datum=WGS84', method = "ngb")
  abunds[is.na(abunds)] <- 0
  abunds_W <- abunds[[2]]
  
  # Calculate relative abundance in each ecoregion
  ecoregions_abund <- terra::extract(abunds_W, ecoregions, sum)
  
  # For each population, get occupied ecoregions 
  points_W_ecoregions <- list()
  pops_name <- vector()
  for(j in 1:length(pops)){
    # Requires at least 3 distinct locations per season
    if(length(which(data_for_analysis2$population == pops[j])) > 0){
      points_W <- SpatialPoints(data_for_analysis2[which(data_for_analysis2$population == pops[j]), 4:5])
      crs(points_W) <-  "+proj=longlat +datum=WGS84"
      points_W_ecoregions[[j]] <- apply(relate(ecoregions, terra::vect(points_W), "contains"), 1, sum)
      pops_name <- c(pops_name, pops[j])
    }
  }
  ecoregions_abund_presence[[i]] <- cbind(ecoregions_abund, do.call(cbind, points_W_ecoregions))
  colnames(ecoregions_abund_presence[[i]]) <- c("ID", "wintering_abund", paste0(pops_name, "_wintering_points"))
}


## Extract worldclim, 2.5 res, 1970-2000 (monthly mean, prec)- https://www.worldclim.org/data/worldclim21.html#
winter <- c("01","02","12")
Temp_files <- list.files("input/distinctiveness/worldclim/wc2.1_2.5m_tavg")
Prec_files <- list.files("input/distinctiveness/worldclim/wc2.1_2.5m_prec")
temp_NB <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[4])) %in% winter], ncol=length(winter))
temp_NB <- terra::rast(apply(temp_NB, 1, function(x) terra::rast(paste0("input/distinctiveness/worldclim/wc2.1_2.5m_tavg/", x)))) # load temperature rasters
temp_NB <- (mean(temp_NB) / 10) - 273.15
prec_NB <- matrix(Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[4])) %in% winter], ncol=length(winter))
prec_NB <- terra::rast(apply(prec_NB, 1, function(x) terra::rast(paste0("input/distinctiveness/worldclim/wc2.1_2.5m_prec/", x)))) # load precipitation rasters
prec_NB <- mean(prec_NB)
crs(temp_NB) <- crs(prec_NB) <- "+proj=longlat +datum=WGS84"


##if you were to use Chelsa climate rasters
## Extract climate for each ecoregion. Climate was downloaded from Chelsa climatologies 1981-2010 (freely available)
winter <- c("01","02","12")
Temp_files <- list.files("~/Dropbox/BGP/Network_modeling_Mignette/LocalAdapt/Manuscript_Revision/distinctiveness/chelsa/tmean")
Prec_files <- list.files("~/Dropbox/BGP/Network_modeling_Mignette/LocalAdapt/Manuscript_Revision/distinctiveness/chelsa/prec")
temp_NB <- matrix(Temp_files[unlist(lapply(strsplit(Temp_files, "_"), function(x) x[3])) %in% winter], ncol=length(winter))
temp_NB <- terra::rast(apply(temp_NB, 1, function(x) terra::rast(paste0("~/Dropbox/BGP/Network_modeling_Mignette/LocalAdapt/Manuscript_Revision/distinctiveness/chelsa/tmean/", x)))) # load temperature rasters
temp_NB <- (mean(temp_NB) / 10) - 273.15
prec_NB <- matrix(Prec_files[unlist(lapply(strsplit(Prec_files, "_"), function(x) x[3])) %in% winter], ncol=length(winter))
prec_NB <- terra::rast(apply(prec_NB, 1, function(x) terra::rast(paste0("~/Dropbox/BGP/Network_modeling_Mignette/LocalAdapt/Manuscript_Revision/distinctiveness/chelsa/prec/", x)))) # load precipitation rasters
prec_NB <- mean(prec_NB)
crs(temp_NB) <- crs(prec_NB) <- "+proj=longlat +datum=WGS84"

# crop around the Americas
temp_NB <- terra::crop(temp_NB, newmap, mask=T)
prec_NB <- terra::crop(prec_NB, newmap, mask=T)

# zscores
temp_mean <- mean(as.vector(temp_NB$mean), na.rm=T)
temp_sd <- sd(as.vector(temp_NB$mean), na.rm=T)
prec_mean <- mean(as.vector(prec_NB$mean), na.rm=T)
prec_sd <- sd(as.vector(prec_NB$mean), na.rm=T)
temp_zscore_NB <- (temp_NB - temp_mean) / temp_sd
prec_zscore_NB <- (prec_NB - prec_mean) / prec_sd

# extract climate for all ecoregions
ecoregions_climate_NB <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions, fun="mean", na.rm=T)
colnames(ecoregions_climate_NB) <- c("ID", "temperature", "precipitation")


## Function to estimate seasonal 2-D climatic niche
nicheDensityRaster <- function(seasonalNiche){
  niche.kernel <- kde2d(seasonalNiche[,1], seasonalNiche[,2], n=50, h=1, lims=c(-2,3, -2,3))
  niche.kernel$z = niche.kernel$z/max(niche.kernel$z)
  niche.raster <- raster(niche.kernel)
  threshold=0; i=0
  while(threshold <= 0.95 * sum(niche.kernel$z)){
    i=i+1
    threshold = threshold + sort(as.vector(niche.raster), decreasing=T)[i]
  }
  niche.raster[which(as.vector(niche.raster) < sort(as.vector(niche.raster), decreasing=T)[i])] = 0
  niche.raster = niche.raster / sum(as.vector(niche.raster))
  return(niche.raster)
}


# Remove the Basin Rockies population of AMRE because of too few wintering samples
ecoregions_abund_presence[[1]] <- ecoregions_abund_presence[[1]][,-7]

# Calculate pairwise niche overlap between wintering populations
climate_distinctiveness_1 <- climate_distinctiveness_2 <- vector()
population_names <- list()
for(k in 1:length(species.names)){
  # Population names
  pops_seas <- do.call(rbind, strsplit(colnames(ecoregions_abund_presence[[k]])[3:ncol(ecoregions_abund_presence[[k]])], "_w"))
  pops_names <- unique(pops_seas[,1])
  population_names[[k]] <- pops_names
  
  # Which ecoregions are occupied by the species
  wintering_ecoregions <- which(ecoregions_abund_presence[[k]]$wintering_abund > 0)
  
  # Percentage of individuals of each population in each ecoregion
  ecoregions_abund_presence_W <- ecoregions_abund_presence[[k]][,3:(2+length(pops_names))]
  ss <- apply(ecoregions_abund_presence_W, 1, sum)
  ss <- ifelse(ss==0, 1, ss)
  ecoregions_abund_presence_W <- apply(ecoregions_abund_presence_W, 2, function(x) x/ss)
  
  # Calculate climate niches (empirical, simulated from ORSIM, simulated from seasonal climate tracking)
  niche_climate_wintering <- list()
  for(j in 3:(2+length(unique(pops_seas[,1])))){
    # which ecoregions are seasonally occupied by a given population
    wintering_pop_ecoregions <- which(ecoregions_abund_presence[[k]][wintering_ecoregions,][,j] > 0)
    
    # Ecoregions weights, based on population relative abundance and the fraction of individuals of the population in each occupied ecoregion 
    ecoregions_weights_wintering <- ecoregions_abund_presence[[k]]$wintering_abund[wintering_ecoregions][wintering_pop_ecoregions] * ecoregions_abund_presence_W[,j-2][wintering_ecoregions][wintering_pop_ecoregions]
    ecoregions_weights_wintering <- ecoregions_weights_wintering / sum(ecoregions_weights_wintering)
    
    # Extract pixels of climate data from climate rasters in each ecoregion
    ecoregions_climate_wintering <- terra::extract(c(temp_zscore_NB, prec_zscore_NB), ecoregions[wintering_ecoregions,][wintering_pop_ecoregions,])
    colnames(ecoregions_climate_wintering) <- c("ID", "temp", "prec") #temp- zcore and prec=zscore
    
    # Resample climate pixels across the species geographical distribution based on ecoregion weights
    ecoregions_climate_resample_wintering <- vector()
    to_sample_wintering <- round(ecoregions_weights_wintering * 10000)
    for(h in 1:length(to_sample_wintering)){
      ecoregions_climate_resample_wintering <- rbind(ecoregions_climate_resample_wintering, ecoregions_climate_wintering[which(ecoregions_climate_wintering$ID == h),][sample(1:length(which(ecoregions_climate_wintering$ID == h)), to_sample_wintering[h], replace=T),])
    }
    ecoregions_climate_resample_wintering <- ecoregions_climate_resample_wintering %>% mutate(population = pops_names[j-2])
    ecoregions_climate_resample_wintering %>% write.table("output/distinctiveness/ecoregions_climate_resample_wintering.all.test4.inclNA.worldclim.oldcomp.txt",row.names=F,quote=F,sep="\t",append=T)
    ecoregions_climate_resample_wintering2<-ecoregions_climate_resample_wintering %>% na.omit()
    # Calculate seasonal climate niches
    niche_climate_wintering[[j-2]] <- nicheDensityRaster(ecoregions_climate_resample_wintering2 %>% dplyr::select(temp, prec))
  }
  #ecoregions_climate_resample_wintering %>% na.omit()
  # Seeasonal climate overlap
  density_wintering <- rasterToPoints(raster::stack(niche_climate_wintering))
  
  distinctiveness_1 <- distinctiveness_2 <- vector()
  for(i in 1:(length(unique(pops_seas[,1]))-1)){
    for(ii in (i+1):length(unique(pops_seas[,1]))){
      # Below are two alternative ways to calculate climate distinctiveness as the opposite of niche overlap. They both give simular results
      # Option 1: using 1 - Schoener's D metric
      distinctiveness_1 <- c(distinctiveness_1, (0.5 * sum(abs(density_wintering[,-c(1,2)][,i] - density_wintering[,-c(1,2)][,ii]))))
      # Option 1: using (1 - overlap) for area with niche density > 0 
      distinctiveness_2 <- c(distinctiveness_2, 1 - length(which(density_wintering[,-c(1,2)][,i] > 0 & density_wintering[,-c(1,2)][,ii] > 0)) / length(which(density_wintering[,-c(1,2)][,i] > 0 | density_wintering[,-c(1,2)][,ii] > 0)))
    }
  }
  #ecoregions_climate_resample_wintering <- rbind(ecoregions_climate_resample_wintering)
  climate_distinctiveness_1[k] <- mean(distinctiveness_1)
  climate_distinctiveness_2[k] <- mean(distinctiveness_2)
  print(k)
}
ecoregions_climate_resample_wintering2 %>% group_by(ID,population) %>% tally()
climate_distinctiveness_1
climate_distinctiveness_2

all<-read.delim("output/distinctiveness/ecoregions_climate_resample_wintering.all.test4.inclNA.worldclim.oldcomp.txt",sep="\t") %>% filter(population!="population") %>% mutate(species=if_else(population=="EST"|population=="PNW_INW"|population=="SW","WIFL",if_else(population=="W_E_Boreal"|population=="PacificCoast"|population=="PNW","SWTH","AMRE"))) %>% na.omit()
all %>% dplyr::group_by(species,population) %>% tally()



ggplot() +
  geom_density(data=all, aes(x = temp,fill = population), alpha = 0.7,adjust = 2) +
 facet_grid(species ~ .,scales='free')+
  labs(
    title = "Temperature Density",
    x = "Temperature (°C)",
    y = "Density"
  ) +
  theme_minimal() + 
  scale_fill_manual(values = c("W_E_Boreal" = "#56B4E9", "NT" = "#56B4E9", "EST" = "#56B4E9","PacificCoast" = "#F0E442","PNW" = "#009E73", "WB" = "#009E73","PNW_INW" = "#009E73","MP" = "#000066","ST" = "#E69F00","SW" = "#E69F00")) + # Assign custom colors
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(panel.border = element_blank(), axis.line = element_line(color = "black")) + 
  theme(aspect.ratio = 1)

ggsave("output/distinctiveness/Temp_zscore.Histogram.Sampling.adjust2.pdf")

ggplot() +
  geom_density(data=all, aes(x = prec,group=population,fill = population), alpha = 0.7,adjust = 2) +
  facet_grid(species ~ .,scales='free')+
  labs(
    title = "Precipitation Density",
    x = "Precipiation (mm)",
    y = "Density"
  ) +
  theme_minimal() + 
  scale_fill_manual(values = c("W_E_Boreal" = "#56B4E9", "NT" = "#56B4E9", "EST" = "#56B4E9","PacificCoast" = "#F0E442","PNW" = "#009E73", "WB" = "#009E73","PNW_INW" = "#009E73","MP" = "#000066","ST" = "#E69F00","SW" = "#E69F00")) + # Assign custom colors
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank(), axis.line = element_line(color = "black")) + theme(aspect.ratio = 1) + xlim(-2,8)

ggsave("output/distinctiveness/Precip_zscore.Histogram.Sampling.adjust2.pdf")


#install.packages("overlapping")
library(overlapping)
swth_t1<-all %>% filter(species=="SWTH") %>% filter(population=="W_E_Boreal") %>% dplyr::select(temp) %>% as.data.frame()
swth_t1n <- as.numeric(swth_t1$temp)
swth_kde_t1 <- density(swth_t1n, kernel = "gaussian",bw=0.25)
range(swth_t1$temp)

swth_t2<-all %>% filter(species=="SWTH") %>% filter(population=="PNW")%>% dplyr::select(temp)%>% as.data.frame()
swth_t2n <- as.numeric(swth_t2$temp)
range(swth_t2$temp)
swth_kde_t2 <- density(swth_t2n, kernel = "gaussian",bw=0.25)

swth_t3<-all %>% filter(species=="SWTH") %>% filter(population=="PacificCoast")%>% dplyr::select(temp)%>% as.data.frame()
swth_t3n <- as.numeric(swth_t3$temp)
swth_kde_t3 <- density(swth_t3n, kernel = "gaussian",bw=0.25)

#Plot the smoothed gaussian kernel density of temp for SWTH
pdf("output/distinctiveness/SWTH.temp_zscore.kde_gaussian.25bw.base_plot.worldclim.pdf")
par(bty = "l")
plot(swth_kde_t2, main = "SWTH Gaussian Kernel Density Estimates of Temp",
     xlab = "Temperature (°C)", ylab = "Density", col ="#009E73" , lwd = 2,type = "l")

# 4. Add the other two densities
lines(swth_kde_t3, col = "#F0E442",lwd = 2)
lines(swth_kde_t1, col = "#56B4E9",lwd = 2)
polygon(swth_kde_t2, col = rgb(0,0.62,0.45, alpha = 0.6))
polygon(swth_kde_t3, col = rgb(0.94,0.89,0.26, alpha = 0.6))
polygon(swth_kde_t1, col = rgb(0.34,0.71,0.91, alpha = 0.6))
dev.off()

##Use code below to get rbg color values for plot above
rgb_values <- col2rgb("#F0E442")
# Print the RGB values
print(rgb_values) #then divide by 255

#create list of numerical distribution to estimate gaussian kernel dist overlap (no smoothing)

x <- list(W_E_Boreal = swth_t1n, PNW = swth_t2n, PacCoast = swth_t3n)
overlap_resultsT <- overlap(x,plot = TRUE)

#overlap with the smoothing
overlap_resultsTs <- ovmult(x,bw=0.25)
ggsave("output/distinctiveness/SWTH.overlap.plot.pdf")

SWTH_distinctT<- 1-overlap_resultsT$OV
SWTH_distinctTs<- 1-overlap_resultsTs$OV

print(paste("SWTH distinct area for temp, bw 0.25:", SWTH_distinctTs))


swth_p1<-all %>% filter(species=="SWTH") %>% filter(population=="W_E_Boreal") %>% dplyr::select(prec) %>% as.data.frame()
dim(swth_p1)
swth_p1n <- as.numeric(swth_p1$prec)
swth_kde_p1 <- density(swth_p1n, kernel = "gaussian",bw=0.25)
range(swth_p1$prec)

swth_p2<-all %>% filter(species=="SWTH") %>% filter(population=="PNW")%>% dplyr::select(prec)%>% as.data.frame()
swth_p2n <- as.numeric(swth_p2$prec)
swth_kde_p2 <- density(swth_p2n, kernel = "gaussian",bw=0.25)
range(swth_p2$prec)

swth_p3<-all %>% filter(species=="SWTH") %>% filter(population=="PacificCoast")%>% dplyr::select(prec)%>% as.data.frame()
swth_p3n <- as.numeric(swth_p3$prec)
swth_kde_p3 <- density(swth_p3n, kernel = "gaussian",bw=0.25)
range(swth_p3$prec)

#This is overlap, but with gaussian distribution, but not smoothing (so no)
p <- list(W_E_Boreal = swth_p1n, PNW = swth_p2n, PacCoast = swth_p3n)
overlap_resultsP <- overlap(p,kernel="gaussian",plot = TRUE)
#smoothing=yes
overlap_resultsPs <- ovmult(p,kernel="gaussian",bw=0.25)

ggsave("output/distinctiveness/SWTH.precip.overlap.plot.no_smoothing.pdf")

SWTH_distinctP<- 1-overlap_resultsP$OV
SWTH_distinctPs<- 1-overlap_resultsPs$OV
SWTH_distinctP
SWTH_distinctPs

print(paste("SWTH distinct area for precip, bw 0.25:", SWTH_distinctPs))

pdf("output/distinctiveness/SWTH.precip_zscore.kde_gaussian.25bw.base_plot.worldclim.NArm_WB.pdf")
par(bty = "l")
plot(swth_kde_p3, main = "SWTH Gaussian Kernel Density Estimates of Precip",
     xlab = "Precipitation (mm)", ylab = "Density", col ="#F0E442" , lwd = 2,type = "l")

# 4. Add the other two densities
lines(swth_kde_p2, col = "#009E73",lwd = 2)
lines(swth_kde_p1, col = "#56B4E9",lwd = 2)
polygon(swth_kde_p3, col = rgb(0.94,0.89,0.26, alpha = 0.6))
polygon(swth_kde_p2, col = rgb(0,0.62,0.45, alpha = 0.6))
polygon(swth_kde_p1, col = rgb(0.34,0.71,0.91, alpha = 0.6))
dev.off()


##WIFL
wifl_t1<-all %>% filter(species=="WIFL") %>% filter(population=="PNW_INW") %>% dplyr::select(temp) %>% as.data.frame()
wifl_t1n <- as.numeric(wifl_t1$temp)
wifl_kde_t1 <- density(wifl_t1n, kernel = "gaussian",bw=0.25)
dim(wifl_t1)
range(wifl_t1$temp)

wifl_t2<-all %>% filter(species=="WIFL") %>% filter(population=="EST")%>% dplyr::select(temp)%>% as.data.frame()
wifl_t2n <- as.numeric(wifl_t2$temp)
wifl_kde_t2 <- density(wifl_t2n, kernel = "gaussian",bw=0.25)
range(wifl_t2$temp)

wifl_t3<-all %>% filter(species=="WIFL") %>% filter(population=="SW")%>% dplyr::select(temp)%>% as.data.frame()
wifl_t3n <- as.numeric(wifl_t3$temp)
wifl_kde_t3 <- density(wifl_t3n, kernel = "gaussian",bw=0.25)
range(wifl_t3$temp)

pdf("output/distinctiveness/WIFL.temp_zscore.kde_gaussian.25bw.base_plot.worldclim.pdf")
par(bty = "l")
plot(wifl_kde_t2, main = "Three Gaussian Kernel Density Estimates of Temp",
     xlab = "Temperature (°C)", ylab = "Density", col ="#56B4E9" , lwd = 2,type = "l") #EST #blue 

# 4. Add the other two densities
lines(wifl_kde_t1, col = "#009E73",lwd = 2) #PNW-INW #009E73 #green
lines(wifl_kde_t3, col = "#E69F00",lwd = 2) #SW orange
polygon(wifl_kde_t2, col = rgb(0.34,0.71,0.91, alpha = 0.6))
polygon(wifl_kde_t1, col = rgb(0,0.62,0.45, alpha = 0.6))
polygon(wifl_kde_t3, col = rgb(0.90,0.62,0, alpha = 0.6))
dev.off()

#Overlap with gaussian kernel density (so ok), no smoothing
x <- list(PNW_INW = wifl_t1n, EST= wifl_t2n, SW = wifl_t3n)
Woverlap_resultsT <- overlap(x,kernel = "gaussian", plot = TRUE)
#with smoothing
Woverlap_resultsTs <- ovmult(x,kernel="gaussian",bw=0.25)

ggsave("output/distinctiveness/WIFL.temp_zscore.overlap.plot.pdf")

WIFL_distinctT<- 1-Woverlap_resultsT$OV
WIFL_distinctTs<- 1-Woverlap_resultsTs$OV

print(paste("WIFL distinct area with smoothing, bw 0.25:", WIFL_distinctTs))

wifl_p1<-all %>% filter(species=="WIFL") %>% filter(population=="PNW_INW") %>% dplyr::select(prec) %>% as.data.frame()
wifl_p1n <- as.numeric(wifl_p1$prec)
wifl_kde_p1 <- density(wifl_p1n, kernel = "gaussian",bw=0.25)
range(wifl_p1$prec)

wifl_p2<-all %>% filter(species=="WIFL") %>% filter(population=="EST")%>% dplyr::select(prec)%>% as.data.frame()
wifl_p2n <- as.numeric(wifl_p2$prec)
wifl_kde_p2 <- density(wifl_p2n, kernel = "gaussian",bw=0.25)
range(wifl_p2$prec)

wifl_p3<-all %>% filter(species=="WIFL") %>% filter(population=="SW")%>% dplyr::select(prec)%>% as.data.frame()
wifl_p3n <- as.numeric(wifl_p3$prec)
wifl_kde_p3 <- density(wifl_p3n, kernel = "gaussian",bw=0.25)
range(swth_p3$prec)

#Yes gaussian kernel distribution (so yes- this is guassian kernel density)
p <- list(PNW_INW = wifl_p1n, EST = wifl_p2n, SW = wifl_p3n)
Woverlap_resultsP <- overlap(p,kernel="gaussian",plot = TRUE)
Woverlap_resultsPs <- ovmult(p,kernel="gaussian",bw=0.25)

ggsave("output/distinctiveness/WIFL.precip_zscore.overlap.plot.pdf")

WIFL_distinctP<- 1-Woverlap_resultsP$OV
WIFL_distinctPs<- 1-Woverlap_resultsPs$OV

print(paste("WIFL distinct area no smooth:", WIFL_distinctP))
print(paste("WIFL distinct area with smoothing:", WIFL_distinctPs)) #use this

pdf("output/distinctiveness/WIFL.precip_zscore.kde_gaussian.25bw.base_plot.worldclim.pdf")
par(bty = "l")
plot(wifl_kde_p3, main = "Three Gaussian Kernel Density Estimates of Precip",
     xlab = "Precipitation (mm)", ylab = "Density", col ="#E69F00" , lwd = 2,type = "l") #SW orange

# 4. Add the other two densities
lines(wifl_kde_p1, col = "#009E73",lwd = 2) #PNW-INW #009E73 #green
lines(wifl_kde_p2, col = "#56B4E9",lwd = 2) #EST #blue
polygon(wifl_kde_p3, col = rgb(0.90,0.62,0, alpha = 0.6))
polygon(wifl_kde_p1, col = rgb(0,0.62,0.45, alpha = 0.6))
polygon(wifl_kde_p2, col = rgb(0.34,0.71,0.91, alpha = 0.6))
dev.off()

##Use code below to get rbg color values for plot above
rgb_values <- col2rgb("#E69F00")
# Print the RGB values
print(rgb_values) #then divide by 255


##AMRE
amre_t1<-all %>% filter(species=="AMRE") %>% filter(population=="MP") %>% dplyr::select(temp) %>% as.data.frame()
amre_t1n <- as.numeric(amre_t1$temp)
amre_kde_t1 <- density(amre_t1n, kernel = "gaussian",bw=0.25)
dim(amre_t1)
range(amre_t1$temp)

amre_t2<-all %>% filter(species=="AMRE") %>% filter(population=="NT")%>% dplyr::select(temp)%>% as.data.frame()
amre_t2n <- as.numeric(amre_t2$temp)
amre_kde_t2 <- density(amre_t2n, kernel = "gaussian",bw=0.25)
range(amre_t2$temp)

amre_t3<-all %>% filter(species=="AMRE") %>% filter(population=="ST")%>% dplyr::select(temp)%>% as.data.frame()
amre_t3n <- as.numeric(amre_t3$temp)
amre_kde_t3 <- density(amre_t3n, kernel = "gaussian",bw=0.25)

amre_t4<-all %>% filter(species=="AMRE") %>% filter(population=="WB")%>% dplyr::select(temp)%>% as.data.frame()
amre_t4n <- as.numeric(amre_t4$temp)
amre_kde_t4 <- density(amre_t4n, kernel = "gaussian",bw=0.25)
range(amre_t4$temp)

 #"NT" = "#56B4E9", 2
 #"WB" = "#009E73", 4
 #"MP" = "#000066", 1
 #"ST" = "#E69F00",  3 # Assign custom colors
  
##Plot in base R (on 1 plot)
pdf("output/distinctiveness/AMRE.temp_zscore.kde_gaussian.25bw.base_plot.worldclim.pdf")
par(bty = "l")
plot(amre_kde_t2, main = "Four Gaussian Kernel Density Estimates of Temp",
     xlab = "Temperature (°C)", ylab = "Density", col ="#56B4E9" , lwd = 2,type = "l") #NT #blue ok

# 4. Add the other two densities
lines(amre_kde_t1, col = "#000066",lwd = 2) #MP #000066 #navy
lines(amre_kde_t4, col = "#009E73",lwd = 2) #WB 009E73 #green
lines(amre_kde_t3, col = "#E69F00",lwd = 2) #ST orange
polygon(amre_kde_t2, col = rgb(0.34,0.71,0.91, alpha = 0.6))
polygon(amre_kde_t1, col = rgb(0,0,0.4, alpha = 0.6))
polygon(amre_kde_t4, col = rgb(0,0.62,0.45, alpha = 0.6))
polygon(amre_kde_t3, col = rgb(0.90,0.62,0, alpha = 0.6))
dev.off()

##Use code below to get rbg color values for plot above
rgb_values <- col2rgb("#009E73")
# Print the RGB values
print(rgb_values) #then divide by 255

atn <- list(MP = amre_t1n, NT= amre_t2n, ST = amre_t3n,WB = amre_t4n)

Aoverlap_resultsT <- overlap(atn,kernel="gaussian",plot = TRUE)
Aoverlap_resultsTs <- ovmult(atn,kernel="gaussian",bw=0.25)

ggsave("output/distinctiveness/AMRE.temp_zscore.overlap.plot.pdf")

amre_distinctT<- 1-Aoverlap_resultsT$OV
amre_distinctTs<- 1-Aoverlap_resultsTs$OV

print(paste("AMRE distinct area no smooth:", amre_distinctT))
print(paste("AMRE distinct area with smoothing:", amre_distinctTs))

#Precip for AMRE
amre_p1<-all %>% filter(species=="AMRE") %>% filter(population=="MP") %>% dplyr::select(prec) %>% as.data.frame()
amre_p1n <- as.numeric(amre_p1$prec)
amre_kde_p1 <- density(amre_p1n, kernel = "gaussian",bw=0.25)
range(amre_p1$prec)

amre_p2<-all %>% filter(species=="AMRE") %>% filter(population=="NT")%>% dplyr::select(prec)%>% as.data.frame()
amre_p2n <- as.numeric(amre_p2$prec)
amre_kde_p2 <- density(amre_p2n, kernel = "gaussian",bw=0.25)
range(amre_p2$prec)

amre_p3<-all %>% filter(species=="AMRE") %>% filter(population=="ST")%>% dplyr::select(prec)%>% as.data.frame()
amre_p3n <- as.numeric(amre_p3$prec)
amre_kde_p3 <- density(amre_p3n, kernel = "gaussian",bw=0.25)
range(amre_p3$prec)

amre_p4<-all %>% filter(species=="AMRE") %>% filter(population=="WB")%>% dplyr::select(prec)%>% as.data.frame()
amre_p4n <- as.numeric(amre_p4$prec)
amre_kde_p4 <- density(amre_p4n, kernel = "gaussian",bw=0.25)
range(amre_p4$prec)

ap <- list(MP = amre_p1n, NT = amre_p2n, ST = amre_p3n, WB= amre_p4n)
Aoverlap_resultsP <- overlap(ap,kernel="gaussian" ,plot = TRUE)
Aoverlap_resultsPs <- ovmult(ap,kernel="gaussian",bw=0.25)

ggsave("output/distinctiveness/AMRE.precip_zscore.overlap.25nbin.plot.pdf")

amre_distinctP<- 1-Aoverlap_resultsP$OV
amre_distinctPs<- 1-Aoverlap_resultsPs$OV

print(paste("AMRE distinct area no smooth:", amre_distinctP))
print(paste("AMRE distinct area with smoothing:", amre_distinctPs))


pdf("output/distinctiveness/AMRE.precip_zscore.kde_gaussian.25bw.base_plot.worldclim.pdf")
par(bty = "l")
plot(amre_kde_p4, main = "Four Gaussian Kernel Density Estimates of Precip",
     xlab = "Precipitation (mm)", ylab = "Density", col ="#009E73" , lwd = 2,type = "l") #NT #blue ok

# 4. Add the other two densities
lines(amre_kde_p1, col = "#000066",lwd = 2) #MP #000066 #navy
lines(amre_kde_p2, col = "#56B4E9",lwd = 2) #WB 009E73 #green
lines(amre_kde_p3, col = "#E69F00",lwd = 2) #ST orange
polygon(amre_kde_p4, col = rgb(0,0.62,0.45, alpha = 0.6))
polygon(amre_kde_p2, col = rgb(0.34,0.71,0.91, alpha = 0.6))
polygon(amre_kde_p1, col = rgb(0,0,0.4, alpha = 0.6))
polygon(amre_kde_p3, col = rgb(0.90,0.62,0, alpha = 0.6))
dev.off()


