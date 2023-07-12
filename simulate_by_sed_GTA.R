# Seed localities across N. America randomly within sediment.
# Do all spatial operations with raster data.
# Run simulations in parallel.
# This script exports the DS3 and DS5 data files that are used 
# in tau and beta analysis.

library(maptools)
library(sp)
library(PBSmapping)
library(rgeos)
library(rgdal)
library(raster)
library(geosphere)
library(vegan)
library(tidyr)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

# save names to put packages on all cores later
pkgs <- c('sp','raster','PBSmapping','vegan','tidyr') 

# set up the parameters of the analysis
n_iter <- 200
grid_res <- 0.5 # 0.5, 1, 5
sites <- round(exp(2:7))
# a natural exponential series

# read in data
spp_brik <- brick(paste0('Data/spatial_data/NAm_mammal_ranges_', grid_res,'degree_any_overlap.grd'))
spp <- names(spp_brik)
sed_brik <- brick(paste0('Data/spatial_data/stage_sediment_', grid_res,'degree_any_overlap.grd'))
bins <- names(sed_brik)
# convert raster brick to list of rasters through which to iterate
sed_l <- lapply(1:length(bins), function(x) sed_brik[[x]])

# sp-level function to calculate range size
# nest this function within the stage-level function, 'sim'
range_sizer <- function(sp, spp_brik, loc_pts){
  sp_map <- spp_brik[[sp]]
  prj <- proj4string(loc_pts)
  if (proj4string(sp_map) != prj){
    stop('wrong projection')
  }
  pres_abs <- sp_map[loc_pts,]
  pres_pts <- loc_pts[pres_abs==1] 
  
  if (length(pres_pts)==0){
    # no fossil sites within species range
    out <- c(sp, rep(NA, 4)) # 5
      
  } else {
    # case where species range overlaps with sampling points:
    # calculate all range metrics
    
    gcell <- length(pres_pts)
    
    ll_prj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    pres_lat <- spTransform(pres_pts,ll_prj)
    lat <- pres_lat@coords[,2]
    latrange <- max(lat) - min(lat)
    
    if (gcell < 2) { 
      convex_area <- gcd <- NA  # mst_length <- 
    } else {
      gcdists <- spDists(pres_pts) 
      gcd <- max(gcdists) / 1000
      
  #    mst_sp <- spantree(gcdists)
  #    mst_length <- sum(mst_sp$dist) / 1000
      
      pres_coords <- pres_pts@coords
      colnames(pres_coords) <- c('X','Y')
      convex <- calcConvexHull(pres_coords)
      convex_area <- calcArea(convex)$area / (10^6)
    }
    
    out <- c(sp, convex_area, gcd, gcell, latrange) # mst_length
  }
return(out)
}

# Stage-level function
sim <- function(sed, s, spp, spp_brik){
  # randomly pick sites to sample within sediment area
  sed_cells <- which(values(sed)==1)
  round(s*length(sed_cells))
  loc_cells <- unique(sample(sed_cells, s, replace = TRUE))
  loc_coords <- xyFromCell(sed, loc_cells)
  sed_prj <- proj4string(sed)
  loc_pts <- SpatialPoints(loc_coords, proj4string = CRS(sed_prj))
  
  # reconstruct range size from simulated fossil sites
  spp_ranges <- sapply(spp, range_sizer, spp_brik=spp_brik, loc_pts=loc_pts)
  spp_df <- data.frame(t(spp_ranges))
  colnames(spp_df) <- c('species','chull','gcd','gcell','latrange') # 'mst'
  out <- gather(spp_df, metric, value, c(chull, gcd, gcell, latrange)) #mst, 
  as.matrix(out)
}

# Iterate and condense output of stage-level function
iter_sim <- function(n_iter, sed, s, spp, spp_brik){
  iter_array <- replicate(n_iter, simplify=FALSE,
                          sim(sed=sed, s=s, spp=spp, spp_brik=spp_brik)
                          )
  
  tmplt <- iter_array[[1]][,1:2]
  iter_ranges <- sapply(iter_array, function(x) x[,3])
  bin <- names(sed)
  out <- cbind(tmplt, stage=bin, sites=s, iter_ranges)
  out <- data.frame(out)
  colnames(out)[5:ncol(out)] <- paste0('rep', 1:n_iter)
  return(out)
}

# Output 'true range size' data for each species and stage

full_map <- spp_brik[[1]]
all_coords <- xyFromCell(full_map, 1:ncell(full_map))
prj <- proj4string(full_map)
all_pts <- SpatialPoints(all_coords, proj4string = CRS(prj))

spp_ranges <- sapply(spp, range_sizer, spp_brik=spp_brik, loc_pts=all_pts)
spp_df <- data.frame(t(spp_ranges))
colnames(spp_df) <- c('species','chull','gcd','gcell','latrange') # 'mst',
# exclude species with too few sites to calculate convex hull
scarce <- which(spp_df$gcell %in% c('1','2'))
spp_df <- spp_df[-scarce,]
spp <- as.character(spp_df$species)
sp_nm <- paste0('Data/DS3_IUCN_range_data_',grid_res,'degree_raster.csv')
write.csv(spp_df, sp_nm, row.names=FALSE)

# apply over all Paleozoic bins and number of sites
# 2.5 hours for 100 replicates

ncores <- detectCores() - 1
pt1 <- proc.time()
registerDoParallel(ncores)

fin <- foreach(sed=sed_l, .packages=pkgs, .combine=rbind, .inorder=FALSE) %:% 
  foreach(s=sites, .packages=pkgs, .combine=rbind, .inorder=FALSE) %dopar% 
  iter_sim(n_iter=n_iter, sed=sed, s=s, spp=spp, spp_brik=spp_brik)

stopImplicitCluster()
pt2 <- proc.time()
(pt2-pt1)/60

fin$sites <- as.numeric(as.character(fin$sites))
fin$sim_id <- paste0(substr(fin$stage, 1, 2), sprintf("%03d", fin$sites), fin$metric)

# save memory: remove rows where all values are NA (i.e. species not sampled in the sim)
rangeCols <- grep('rep', colnames(fin))
sumRange <- rowSums(fin[,rangeCols], na.rm = TRUE) 
fin <- fin[sumRange != 0,]

fin_nm <- paste0('Data/DS5_simulation_range_data_',grid_res,'degree_by_sed.csv')
write.csv(fin, fin_nm, row.names = FALSE)
