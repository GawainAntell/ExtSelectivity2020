library(beepr)
library(stringr)
library(raster)
library(sp)
library(rgdal)

grid_res <- 0.5 

# Equal-area projected template for rasterizing polygons
us_bounds <- c(-135, -65, 22, 50)
prj <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +
  ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# "+proj=eck6 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
r_ll <- raster(resolution=grid_res, ext=extent(us_bounds))
r_prj <- projectRaster(r_ll, crs=prj)
values(r_prj) <- 1:ncell(r_prj)

# North American mammal ranges

# Remove species not found within continental US extent
mam_data <- readOGR("Data/shapefiles/TERRESTRIAL_MAMMALS.shp")
mam_clipped <- crop(mam_data, extent(us_bounds))

# use exactly the same 342 species as by-species approach
# e.g. exclude Ursus americanus
vetted <- read.csv('Data/DS2_IUCN_range_data_vector.csv',stringsAsFactors=FALSE)
spp <- gsub('_',' ', vetted$species)

# Format N. American mammal polygons as a brick:
spp_brik <- brick()
for (sp_nm in spp) {  
  sp_poly <- which(mam_clipped$scientific == sp_nm)
  sp_spdf <- mam_clipped[sp_poly,]
  sp_spdf <- spTransform(sp_spdf, prj)
  sp_cells_l <- raster::extract(r_prj, sp_spdf, weights=TRUE, small=TRUE)
  sp_cells <- unique(unlist(sp_cells_l))
  sp_r <- r_prj
  values(sp_r) <- 0
  values(sp_r)[sp_cells] <- 1
  names(sp_r) <- sp_nm
  spp_brik <- brick(list(spp_brik, sp_r))
} 
beep('treasure')
# native 'raster' format will preserve layers names when read back in to R
# but also export a header file in order to open file in other app, e.g. QGIS
sp_nm <- paste0('Data/NAm_mammal_ranges_',grid_res,'degree_any_overlap')
writeRaster(spp_brik, sp_nm, format='raster', overwrite=TRUE)
hdr(spp_brik, filename=sp_nm, format='ENVI')

# Format Paleozoic sediment polygons as a brick:
all_sed_fls <- list.files('Data/Stage_shp_files', recursive=TRUE)
shp_fl_pos <- grep('.shp', all_sed_fls)
xml_fl_pos <- grep('.shp.xml', all_sed_fls)
shp_fl_pos <- setdiff(shp_fl_pos, xml_fl_pos)
shp_fls <- all_sed_fls[shp_fl_pos]
shp_fls_dsn <- paste0('Data/Stage_shp_files/', shp_fls)
bins <- sapply(shp_fls, function(x) str_split(x, '_')[[1]][2] )
n_bins <- length(bins)
sed_brik <- brick()
for (j in 1:n_bins){
  sed_j <- readOGR(shp_fls_dsn[j], verbose = FALSE)
  sed_prj <- spTransform(sed_j, prj)
  sed_cells_l <- raster::extract(r_prj, sed_prj, weights=TRUE, small=TRUE)
  sed_cells <- unique(unlist(sed_cells_l))
  sed_r <- r_prj
  values(sed_r) <- 0
  values(sed_r)[sed_cells] <- 1
  names(sed_r) <- bins[j]
  sed_brik <- brick(list(sed_brik, sed_r))
}
sed_nm <- paste0('Data/stage_sediment_',grid_res,'degree_any_overlap')
writeRaster(sed_brik, sed_nm, format='raster', overwrite=TRUE)
hdr(sed_brik, filename=sed_nm, format='ENVI')
