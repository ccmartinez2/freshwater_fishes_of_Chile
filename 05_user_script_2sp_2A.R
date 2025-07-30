# -------------------------------------------------------------------------
## Written by Olivier Broennimann 👨‍🔬. Departement of Ecology and Evolution (DEE), University of Lausanne, Switzerland.
#  Modified by: Dr. Jan Arif 🧑‍💻-  OUS Oregon State university Date: 2023
#  Modified by: Cristian Martinez 🧑‍💻-  Department of ecology, UC - Chile. Date: 2024-09-23
#  

# ---------------------------------------------------------------------------------------------------
          ##### Script for niche overlap analysis for Rainbow Trout introduced in Chile🌈 🐠 ####

## free memory 🧠
rm(list = ls())

# -------------------------------------------------------------------------
## Directory 📙 

setwd("C:/freshwater_variables/models/wld_intro_spp/om_na_ch")
getwd()
dir()

# -------------------------------------------------------------------------
#Start
t1=Sys.time()

# -------------------------------------------------------------------------
#Packages 📂 

require(doParallel)
require(tidyverse)
require(devtools)
require(raster)
require(ade4)
require(rgdal)
require(ecospat)
require(rworldmap)
require(sf)
require(terra)
require(pbapply)
require(future.apply)
require(factoextra)
require(RColorBrewer)

# -------------------------------------------------------------------------

                              ###########################
                              #### data preparation #####
                              ###########################

# -------------------------------------------------------------------------
# load environmental data  # get climate data for the two area ☔ 🌐
library(raster)
library(sf)
library(adehabitatHR) 
library(vegan) 

# Load environmental data 
list <- list.files(path = "C:/freshwater_variables/models/wld_intro_spp/om_na_ch", pattern = ".tif$", full.names = TRUE)
subset_list <- list[grep("lc", list, invert = TRUE)] # Subset sin land cover 
clim <- raster::stack(subset_list)

# Convert raster stack to data frame for processing
clim_values <- as.data.frame(raster::getValues(clim))

# Remove outliers
#  outliers limits  (e.x using M&M  of 1.5 * IQR)
outlier_limits <- apply(clim_values, 2, function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  c(lower = q1 - 1.5 * iqr, upper = q3 + 1.5 * iqr)
})

# Delete outliers
clim_values_filtered <- clim_values  # Copy filter
for (i in seq_along(outlier_limits[1, ])) {
  clim_values_filtered[!is.na(clim_values_filtered[, i]) & 
                         (clim_values_filtered[, i] < outlier_limits[1, i] | 
                            clim_values_filtered[, i] > outlier_limits[2, i]), i] <- NA
}

# Convert back to raster stack, keeping ext and res
clim_filtered <- stack()
for (i in 1:ncol(clim_values_filtered)) {
  clim_filtered <- stack(clim_filtered, setValues(clim[[i]], clim_values_filtered[, i]))
}

# Reasigne the env clim obj
clim <- clim_filtered

# Create mask for native area background 
bkg.nam_sf <- st_read("C:/freshwater_variables/models/wld_intro_spp/om_na_ch/nam2.shp")
bkg.chl_sf <- st_read("C:/freshwater_variables/models/wld_intro_spp/om_na_ch/chl2.shp")

# Convert to SpatialPolygons
bkg.nam <- as(bkg.nam_sf, "Spatial")
bkg.chl <- as(bkg.chl_sf, "Spatial")

# Extract environmental data from the rasters for native and introduced areas
clim.bkg.nam <- mask(crop(clim, bbox(bkg.nam)), bkg.nam)
clim.bkg.chl <- mask(crop(clim, bbox(bkg.chl)), bkg.chl)

# -------------------------------------------------------------------------
# load occurrence data  🐠🐠🐠🐠
occ.all<- read_csv("occ.csv")
occ.all

# Convert tibble to data frame
occ.all_df <- as.data.frame(occ.all)

# remove occurrences closer than a minimum distance to each other (remove aggregation)
occ.sp<-ecospat.occ.desaggregation(occ.all,min.dist = 1/60*10)
occ.sp

# Specify which columns contain the coordinates (x and y)
coordinates(occ.all_df) <- ~ x + y
proj4string(occ.all_df) <- crs(clim)

# Extract only the x and y columns (longitude and latitude)
occ_coords <- occ.sp[, c("x", "y")]

#Extract function with the subsetted coordinates
data_occ_all <- data.frame(raster::extract(clim, occ_coords))
data_occ_all<- cbind(data_occ_all,occ.sp)
data_occ_all<- na.omit(data_occ_all)
data_occ_all

### subset env data for species and areas### 
occ.nam<-base::subset(data_occ_all,region=="NAM",select = c(x,y))
occ.chl<-base::subset(data_occ_all,region=="CHL",select = c(x,y))

# -------------------------------------------------------------------------
# Configure parallel 
plan(multisession, workers = 4)

# create sp occurrence dataset by extracting climate variables from the rasters
env_bkg_nam <-na.exclude(data.frame(raster::extract(clim, bkg.nam)))  ###############clim.bkg.nam
env_bkg_chl <-na.exclude(data.frame(raster::extract(clim, bkg.chl)))

env_occ_nam <- na.exclude(data.frame(raster::extract(clim, occ.nam)))
env_occ_chl <- na.exclude(data.frame(raster::extract(clim, occ.chl)))

# -------------------------------------------------------------------------
#check spatial autocorrelation
jpeg("correlograms.jpg", width = 3000, height = 1500, res = 300)

par(mfrow = c(1, 2), mar = c(5, 5, 5, 2))  # Ajustar márgenes: c(inf, izq, sup, der)

#"Correlogram in Native area"
ecospat.mantel.correlogram(
                          dfvar = cbind(occ.nam, env_occ_nam),
                          colxy = 1:2,
                          n = 200,
                          nclass = 10,
                          nperm = 20
)
title("Correlogram in Native area")

# Correlogram in introduced area in Chile
ecospat.mantel.correlogram(
                            dfvar = cbind(occ.chl, env_occ_chl),
                            colxy = 1:2,
                            n = 200,
                            nclass = 10,
                            nperm = 20
)
title("Correlogram in introduced area - Chile")

# Restablish  graph disposition
par(mfrow = c(1, 1))
dev.off()

-----------------------------------------------------------------------
################################
#### niche quantification ###### ✖️➕➖➗🟰
################################

#cool names#
s
env_bkg_nam1<- env_bkg_nam
colnames(env_bkg_nam1) <- names(env_bkg_nam)                               
env_bkg_chl1<- env_bkg_chl
colnames(env_bkg_chl1) <- names(env_bkg_chl)

# -------------------------------------------------------------------------
#calibration of PCA-env 💹

pca.env <-dudi.pca(rbind(env_bkg_nam1,env_bkg_chl1), 
                   center = T, 
                   scale = T, 
                   scannf = F, 
                   nf = 2)

pca_env_nam_chl <-dudi.pca(rbind(env_bkg_nam1,env_bkg_chl1), 
                           center = T, 
                           scale = T, 
                           scannf = F, 
                           nf = 2)

pca_env_chl <-dudi.pca(rbind(env_bkg_chl), 
                       center = T, 
                       scale = T, 
                       scannf = F, 
                       nf = 2)

# -------------------------------------------------------------------------
#################### predict the scores on the PCA axes ####################
# predict the scores on the PCA axes

scores.bkg<- pca.env$li

scores.bkg.nam<- suprow(pca.env,env_bkg_nam1)$lisup
scores.bkg.chl<- suprow(pca.env,env_bkg_chl1)$lisup

scores.occ.nam<- suprow(pca.env,env_occ_nam)$lisup
scores.occ.chl<- suprow(pca.env,env_occ_chl)$lisup

#################### env niches  ####################

# delete  NA values
scores.bkg <- na.omit(scores.bkg)
scores.bkg.nam <- na.omit(scores.bkg.nam)
scores.occ.nam <- na.omit(scores.occ.nam)

scores.occ.chl <- na.omit(scores.occ.chl)
scores.occ.chl <- na.omit(scores.occ.chl)

# Create a .jpg image file with high resolution
jpeg("niche_plots4.jpg", width = 3000, height = 1500, res = 300)

# Set the layout for two plots in one row
par(mfrow = c(1, 2), mar = c(5, 5, 5, 2))  # Adjust margins

#################### Calculate the niche in the native area ####################
zrt_native <- ecospat.grid.clim.dyn(
                                      scores.bkg, scores.bkg.nam,
                                      scores.occ.nam, R = 1000,
                                      kernel.method = "adehabitat",
                                      extend.extent = c(-5, 5, -3, 3))

#extend.extent = c(x_min_extension, x_max_extension, y_min_extension, y_max_extension

# Plot the niche in the native area
ecospat.plot.niche(
                    zrt_native, 
                    title = "RT Native", 
                    name.axis1 = "PC1", 
                    name.axis2 = "PC2", 
                    cor = TRUE)

#################### Calculate the niche in the introduced area (Chile) ####################
zrt_chl <- ecospat.grid.clim.dyn(
                                  scores.bkg, scores.bkg.chl,
                                  scores.occ.chl, R = 1000,
                                  kernel.method = "adehabitat",
                                  extend.extent = c(-5, 5, -3, 3))

# Plot the niche in the introduced area (Chile) with a custom title
ecospat.plot.niche(
                    zrt_chl, 
                    title = "RT Chile", 
                    name.axis1 = "PC1", 
                    name.axis2 = "PC2", 
                    cor = TRUE)

# Reset the graphical layout to the default
par(mfrow = c(1, 1))

# Close the image file to save the plots
dev.off()

#################### Calculate and display the niche overlap between the two areas ####################
niche_overlap <- ecospat.niche.overlap(zrt_native, zrt_chl, cor = TRUE)$D
print(paste("Niche Overlap (D):", niche_overlap))

#################### test of niche equivalency ####################
?ecospat.niche.equivalency.test

# Perform the niche equivalency test
equ <- ecospat.niche.equivalency.test(
                                        zrt_native,
                                        zrt_chl,
                                        rep = 1000,
                                        ncores = 10
                                      )

# Create a .jpg image file with high resolution
jpeg("niche_equivalency_test.jpg", width = 3000, height = 1500, res = 300)

# Plot the overlap test results with a custom title
ecospat.plot.overlap.test(equ, "D", "RT Equivalency Test")

# Close the image file to save the plot
dev.off()

# Print the equivalency test results to the console
print(equ)

#################### test of niche similarity ##################
?ecospat.niche.similarity.test

# Conduct the similarity test for both areas
sim_rt1 <- ecospat.niche.similarity.test(
                                          zrt_native,
                                          zrt_chl,
                                          rep = 1000,
                                          overlap.alternative = "higher",
                                          expansion.alternative = "lower", 
                                          stability.alternative = "higher", 
                                          unfilling.alternative = "lower",
                                          rand.type = 1,  # niches randomly shifted in both areas
                                          ncores = 10
                                        )
print(sim_rt1)  # Print the results

# Conduct the similarity test only in the invaded area
sim_rt2 <- ecospat.niche.similarity.test(
                                          zrt_native,
                                          zrt_chl,
                                          rep = 1000,
                                          overlap.alternative = "higher",
                                          expansion.alternative = "lower", 
                                          stability.alternative = "higher", 
                                          unfilling.alternative = "lower",
                                          rand.type = 2,  # niche randomly shifted only in the invaded area
                                          ncores = 1
                                        )
print(sim_rt2)  # Print the results

# Create a .jpg image file with high resolution
jpeg("niche_similarity_tests.jpg", width = 3000, height = 1500, res = 300)

# Set up the plotting area for two plots in one row
par(mfrow = c(1, 2), mar = c(5, 5, 5, 2))  # Adjust margins

# Plot the similarity test results for both areas
ecospat.plot.overlap.test(sim_rt1, "D", "RT Similarity Test Both Areas")

# Plot the similarity test results only in the invaded area
ecospat.plot.overlap.test(sim_rt2, "D", "RT Similarity Test Only in Invaded Area")

# Close the image file to save the plots
dev.off()

################### test of niche divergence  ##################

# Conduct the niche divergence test
sim_rt_div <- ecospat.niche.similarity.test(
                                              zrt_native,
                                              zrt_chl,
                                              rep = 1000,
                                              intersection = 0,
                                              overlap.alternative = "lower",
                                              rand.type = 2,  # niche randomly shifted only in the invaded area
                                              expansion.alternative = "higher",
                                              stability.alternative = "lower",
                                              unfilling.alternative = "higher",
                                              ncores = 1  # Set to 1 core for this operation
                                            )

# Print the results of the divergence test
print(sim_rt_div)

# Create a .jpg image file with high resolution for the plot
jpeg("niche_divergence_test.jpg", width = 3000, height = 1500, res = 300)

# Plot the results of the divergence test with a custom title
ecospat.plot.overlap.test(sim_rt_div, "D", "RT Divergence Test")

# Close the image file to save the plot
dev.off()

################### plot niche dynamics  ##################

clim1 <- as.matrix(scores.occ.nam[, 1:2]) 
clim2 <- as.matrix(scores.occ.chl[, 1:2])

# Create a .jpg image file with high resolution for the niche plot
jpeg("niche_dynamics_rainbow_trout_3.jpg", width = 3000, height = 1500, res = 500)

# Choose a color palette
colors <- brewer.pal(3, "Set1")  # You can choose other palettes like "Dark2", "Accent", etc.

# Plot the dynamic niche of Rainbow Trout 🖼️🎴
ecospat.plot.niche.dyn(zrt_native, zrt_chl,
                       intersection = 0, 
                       title = "", 
                       name.axis1 = "PC1", 
                       name.axis2 = "PC2", 
                       interest = 1, 
                       quant = 0.0, 
                       col.unf = "lightblue",
                       col.exp = "pink", 
                       col.stab = "gray", 
                       colZ1 = "#009E73", 
                       colZ2 = "#661100",
                       transparency = 50, 
                       xlim = c(-5, 5),     # Adjust x-axis limits
                       ylim = c(-5, 5)      # Adjust y-axis limits
                        )  
# Niche centroids shift
ecospat.shift.centroids(scores.occ.nam, 
                        scores.occ.chl, 
                        clim1 = clim1,
                        clim2 = clim2,
                        col = "black") 



# Close the image file to save the plot
dev.off()


#-------------------------------------------------------------------------------------------------------------------------


# Calculate the niche dynamic index
niche_dyn_index <- ecospat.niche.dyn.index(
                                            zrt_native,
                                            zrt_chl,
                                            intersection = 0
                                          )

# Calculate overlap corrected by availability of background conditions
corrected_overlap <- ecospat.niche.overlap(
                                            zrt_native,
                                            zrt_chl,
                                            cor = TRUE
                                          )

# Calculate uncorrected overlap
uncorrected_overlap <- ecospat.niche.overlap(
                                              zrt_native,
                                              zrt_chl,
                                              cor = FALSE
                                            )

# Optionally print the results of the niche dynamics index and overlaps
print(niche_dyn_index)
print(corrected_overlap)
print(uncorrected_overlap)

##### contribution of original variables ######
# Load necessary libraries
library(factoextra)  # For PCA plotting
library(gridExtra)   # For arranging multiple ggplots
library(ggplot2)     # Ensure ggplot2 is loaded for plotting

# Create PCA contribution plots separately
plot1 <- fviz_pca_var(pca_env_nam_chl,
                      col.var = "contrib",
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE) +
  ggtitle("PCA Contributions in Native Area") +
  theme_minimal()  # Optional: change the theme if needed

plot2 <- fviz_pca_var(pca_env_chl,
                      col.var = "contrib",
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE) +
  ggtitle("PCA Contributions in Introduced Area") +
  theme_minimal()  # Optional: change the theme if needed

# Create a .jpg image file
jpeg("PCA_contributions.jpg", width = 3000, height = 1500, res = 300)

# Combine the plots into a single plot layout
grid.arrange(plot1, plot2, nrow = 1)

# Close the JPEG device
dev.off()

print(equ) # Print niche equivalency test (see Warren et al 2008) based on two species occurrence density grids.
print(sim_rt2)  # Print the  similarity niche test results  niche similarity test (see Warren et al 2008) based on two species occurrence density grids.
print(sim_rt_div) # Print the results of the divergence test

print(niche_dyn_index)  # Print the niche categories of niche dynamics between two species densities 
print(corrected_overlap) # Print niche overlap results
print(uncorrected_overlap)  # Print uncorrected niche overlap results
 

# Results data frame
results <- data.frame(
  Test = c("Similarity Test", "Divergence Test", "Equivalency Test", "Niche Dynamics", "Corrected Overlap", "Uncorrected Overlap"),
  
  # Similarity Test
  D_sim = c(sim_rt2$obs$D, NA, NA, NA, NA, NA),
  I_sim = c(sim_rt2$obs$I, NA, NA, NA, NA, NA),
  Expansion_sim = c(sim_rt2$obs$expansion, NA, NA, NA, NA, NA),
  Stability_sim = c(sim_rt2$obs$stability, NA, NA, NA, NA, NA),
  Unfilling_sim = c(sim_rt2$obs$unfilling, NA, NA, NA, NA, NA),
  p_D_sim = c(sim_rt2$p.D, NA, NA, NA, NA, NA),
  p_I_sim = c(sim_rt2$p.I, NA, NA, NA, NA, NA),
  p_expansion_sim = c(sim_rt2$p.expansion, NA, NA, NA, NA, NA),
  p_stability_sim = c(sim_rt2$p.stability, NA, NA, NA, NA, NA),
  p_unfilling_sim = c(sim_rt2$p.unfilling, NA, NA, NA, NA, NA),
  
  # Divergence Test
  D_div = c(NA, sim_rt_div$obs$D, NA, NA, NA, NA),
  I_div = c(NA, sim_rt_div$obs$I, NA, NA, NA, NA),
  Expansion_div = c(NA, sim_rt_div$obs$expansion, NA, NA, NA, NA),
  Stability_div = c(NA, sim_rt_div$obs$stability, NA, NA, NA, NA),
  Unfilling_div = c(NA, sim_rt_div$obs$unfilling, NA, NA, NA, NA),
  p_D_div = c(NA, sim_rt_div$p.D, NA, NA, NA, NA),
  p_I_div = c(NA, sim_rt_div$p.I, NA, NA, NA, NA),
  p_expansion_div = c(NA, sim_rt_div$p.expansion, NA, NA, NA, NA),
  p_stability_div = c(NA, sim_rt_div$p.stability, NA, NA, NA, NA),
  p_unfilling_div = c(NA, sim_rt_div$p.unfilling, NA, NA, NA, NA),
  
  # Equivalency Test
  D_equ = c(NA, NA, equ$obs$D, NA, NA, NA),
  I_equ = c(NA, NA, equ$obs$I, NA, NA, NA),
  Expansion_equ = c(NA, NA, equ$obs$expansion, NA, NA, NA),
  Stability_equ = c(NA, NA, equ$obs$stability, NA, NA, NA),
  Unfilling_equ = c(NA, NA, equ$obs$unfilling, NA, NA, NA),
  p_D_equ = c(NA, NA, equ$p.D, NA, NA, NA),
  p_I_equ = c(NA, NA, equ$p.I, NA, NA, NA),
  p_expansion_equ = c(NA, NA, equ$p.expansion, NA, NA, NA),
  p_stability_equ = c(NA, NA, equ$p.stability, NA, NA, NA),
  p_unfilling_equ = c(NA, NA, equ$p.unfilling, NA, NA, NA),
  
  # Niche Dynamics
  Expansion_dyn = c(NA, NA, NA, niche_dyn_index$dynamic.index.w["expansion"], NA, NA),
  Stability_dyn = c(NA, NA, NA, niche_dyn_index$dynamic.index.w["stability"], NA, NA),
  Unfilling_dyn = c(NA, NA, NA, niche_dyn_index$dynamic.index.w["unfilling"], NA, NA),
  
  # Corrected Overlap
  D_corrected = c(NA, NA, NA, NA, corrected_overlap$D, NA),
  I_corrected = c(NA, NA, NA, NA, corrected_overlap$I, NA),
  
  # Uncorrected Overlap
  D_uncorrected = c(NA, NA, NA, NA, NA, uncorrected_overlap$D),
  I_uncorrected = c(NA, NA, NA, NA, NA, uncorrected_overlap$I)
)

# Guardar en un archivo .csv
write.csv(results, file = "niche_tests_results.csv", row.names = FALSE)


# Check the working directory to find your file
getwd()  # This will show where the file is saved


