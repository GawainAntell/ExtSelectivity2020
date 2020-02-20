### SAFD code for simulating extant mammal range sizes based on the current distribution of surface fossiliferous sediments
### Only things that should need changing are:

# the working directory (i.e. where you want the results .csvs to be sent)
# the file paths for the US map, IUCN species, and geological .shp files (immediately below the working directory path)
# the parameters of the analysis, i.e. the geological periods, no. of sites to be simulated, no. of iterations, and the size (in deg.) of cells for the grid cell method (clearly labeled with 'KEY PARAMETERS LOCATED HERE')

### Code is currently set up to print the name of the stage currently being run, and a way of tracking progress
### Note that the code also outputs a .csv file 'time_res.csv' detailing the time elapsed (in hours) for each geological stage to be completed

#####################################################################################################################
#####################################################################################################################

# load up the requisite packages

library(maptools)
library(sp)
library(PBSmapping)
library(rgeos)
library(spdep)
library(rgdal)
library(splancs)
library(alphahull)
library(pastecs)
library(plotrix)
library(raster)
library(geosphere)
library(ggplot2)
library(tictoc)

#####################################################################################################################
#####################################################################################################################

# set up a working directory - remember to change file path

setwd("~/Desktop/Projects/Range Size Simulations_w.Geography/experiment_files")

#####################################################################################################################
#####################################################################################################################

# load up the requisite .shp files - remember to change file paths

### US .shp file
US_map<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/Range Size Simulations_w.Geography/US_mapfiles/cb_2014_us_nation_5m.shp", layer = "cb_2014_us_nation_5m")
US_map_rangeproj<-spTransform(US_map, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

### geo period data
# Cambrian<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Cambrian/US_Cambrian_final.shp", layer = "US_Cambrian_final")
# Ordovician<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Ordovician/us_ordovician_final.shp", layer = "us_ordovician_final")
# Silurian<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Silurian/us_silurian_final.shp", layer = "us_silurian_final")
# Devonian<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Devonian/us_devonian_final.shp", layer = "us_devonian_final")
# Carboniferous<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Carboniferous/us_carboniferous.shp", layer = "us_carboniferous")
# Permian<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Permian/us_permian_final.shp", layer = "us_permian_final")
# Triassic<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Triassic/us_triassic_final.shp", layer = "us_triassic_final")
# Jurassic<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Jurassic/us_jurassic_final.shp", layer = "us_jurassic_final")
Cretaceous<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/Range Size Simulations_w.Geography/Contiguous US/Cretaceous/Cretaceous_New_diss.shp", layer = "Cretaceous_New_diss")
# Paleogene<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Paleogene/us_paleogene_final.shp", layer = "us_paleogene_final")
# Neogene<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Neogene/us_neogene_final.shp", layer = "us_neogene_final")
# Quarternary<-readOGR(dsn = "/Users/darrocsa/Desktop/Projects/RangeSizeSimulations_Geography/Contiguous US/Quaternary/us_quaternary_final.shp", layer = "us_quaternary_final")

# get mammal data
mam_data <- readShapeSpatial("/Users/darrocsa/Desktop/Projects/Range Size Simulations_w.Geography/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp")
mam_unique <- unique(mam_data@data$scientific)

#####################################################################################################################
#####################################################################################################################

# load up the IUCN list of N. American mammals/amphibians (N.American_mammals_lim) and extract a vector of names:

mams <- read.csv("~/Desktop/Projects/Range Size Simulations_w.Geography/N.American_mammals_lim.csv")
mam.sp.names <- as.vector(paste(mams$Genus, mams$Species, sep=" "))

# set up an empty list to put subsetted N. American mammals in:
NAtaxa<-list()
counter=1

for (t in 1:length(mam_unique)) {    
    ts<-mam_unique[t]
    tmp<-mam_data[which(mam_data$scientific == ts) , ]
    name<-as.character(tmp$scientific[1])
    # print(name)   
    if (name %in% mam.sp.names == TRUE){    	
    	NAtaxa[[counter]]<-tmp
    	counter=counter+1    	
    }  
    # print(counter)    
}    

### see how many species we have in there...
# length(NAtaxa) # 377
all_species<-NAtaxa
l.species<-length(all_species)

# get a vector species names
sp.names<-c()

for (smut in 1:length(NAtaxa)){
	t.smut<-NAtaxa[[smut]]
	nm<-t.smut$scientific[1]
	sp.names<-c(sp.names, as.character(nm))
}

######################################## KEY PARAMETERS LOCATED HERE ################################################
#####################################################################################################################

# set up the parameters of the analysis

# define no. of periods and species
all_periods<-list(Cretaceous)
p.names <- c("Cretaceous")
# all_periods<-list(Cambrian, Ordovician, Silurian, Devonian, Carboniferous, Permian, Triassic, Jurassic, Cretaceous, Paleogene, Neogene, Quarternary)
# p.names <- c("Cambrian", "Ordovician", "Silurian", "Devonian", "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quarternary")

# define iterations
iterations<-c(100)

# define vector for no. of sites simulated
sites<-c(200)

# define size of sampling grid for grid cell method
grid_res<-5

#####################################################################################################################
#####################################################################################################################

# miscellaneous things that should be automatically sorted out...

# define lengths
l.periods<-length(all_periods)
l.species<-length(all_species)
l.sites<-length(sites)

# make a raster for the grid cell method:
x.raster<-raster(ncol=360/grid_res, nrow=180/grid_res, xmn=-180, xmx=180, ymn=-90, ymx=90)
x.raster<-setValues(x.raster, c(1:ncell(x.raster)))

# define lengths
l.periods<-length(all_periods)
l.species<-length(all_species)
l.sites<-length(sites)

#####################################################################################################################
#####################################################################################################################

# set up a quick results table for recording the time elapsed taken to complete each stage
time_res<-array(NA, dim=c(length(all_periods), l.sites))
colnames(time_res)<-c(sites)
rownames(time_res)<-c(p.names)

#####################################################################################################################
#####################################################################################################################

# Analysis

for (i in 1:l.periods){
	
	per<-all_periods[[i]]
	c.period<-spTransform(per, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
	
	# identify the period and give ourselves a check
	p.name<-p.names[i]
	# print(p.name)
	
		for (s in 1:l.sites){
		
		# start the clock
		a<-tic()
		
		# establish no. sites	
		t.sites<-sites[s]
					
		### make some results arrays
		# make a results table for lat. range
		res_array_latrange<-array(NA, dim=c(iterations, l.species))
		colnames(res_array_latrange)<-c(as.character(sp.names))
		lr.s.name<-paste(p.name, "latrange", t.sites, "sites", "csv", sep=".")
	
		# make a results table for GCD
		res_array_gcd<-array(NA, dim=c(iterations, l.species))
		colnames(res_array_gcd)<-c(as.character(sp.names))
		gcd.s.name<-paste(p.name, "gcd", t.sites, "sites", "csv", sep=".")
		
		# make a results table for convex hull
		res_array_chull<-array(NA, dim=c(iterations, l.species))
		colnames(res_array_chull)<-c(as.character(sp.names))
		chull.s.name<-paste(p.name, "chull", t.sites, "sites", "csv", sep=".")
		
		# make a results table for grid cell
		res_array_gcell<-array(NA, dim=c(iterations, l.species))
		colnames(res_array_gcell)<-c(as.character(sp.names))
		gcell.s.name<-paste(p.name, "gcell", t.sites, "sites", grid_res, "deg", "csv", sep=".")
		
		# make a results table for grid cell percentage
		# res_array_gcell_perc <-array(NA, dim=c(iterations, l.species))
		# colnames(res_array_gcell_perc)<-c(as.character(sp.names))
		# gcell.s.perc.name<-paste(p.name, "gcell", t.sites, "sites", grid_res, "percentage", "deg", "csv", sep=".")
	
			for (j in 1:l.species){
		
			c.species<-all_species[[j]]
			projection(c.species)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"	# get the range in the right projection
				
			# cut the polygons
			species_cropped <-gIntersection(c.period, c.species)
		
			sp.name<-sp.names[j]
			# print(sp.name)
			
			# print a check string
			check<-paste(p.name, t.sites, sp.name, sep="    ")
			print(check)
		
			if (is.null(species_cropped)) { 			
				
				res_array_latrange[,j] <- c(NA)
				res_array_gcd[,j] <- c(NA)
				res_array_chull[,j] <- c(NA)
				res_array_gcell[,j] <- c(NA)
			
			} else {		
	
				for (k in 1:iterations){
			
				# print(k)
				
				# simulate a fossil record
				extinction_record<-spsample(species_cropped, sites, type="random", iter=1000000)
			
				### start with lat. range
				latlong_record<-cbind(extinction_record$x, extinction_record$y)
				max.lat<-max(latlong_record[,2])
				min.lat<-min(latlong_record[,2])
				sim.lat.range<-abs(max.lat - min.lat)
			
				# save the result
				res_array_latrange[k,j]<-sim.lat.range
				
				### then GCD
				# set up an empty vector and organize the coordinates
				gcdists<-vector()
				gc_record<-cbind(extinction_record$x, extinction_record$y)
		
				# then create a nested loop to perform a great circle cal. (spDistsN1) on each point, and tape the results to the vector
				for (row in 1:nrow(gc_record)){
					now.row<-gc_record[row,]
					gcs<-spDistsN1(gc_record, now.row, longlat = TRUE)
					gcdists<-c(gcdists, gcs)						
				}

				# now need to get the maximum of all these possible combinations, and work out how close it is to the actual maximum great circle distance of the range
				max.gcdist<-max(gcdists)
				res_array_gcd[k,j]<-max.gcdist
				
				### then convex hull
				convex_record_alt1<-cbind(extinction_record$x, extinction_record$y)
				convex_record_data1<-as.data.frame(convex_record_alt1)
				convex_spatial1<-SpatialPointsDataFrame(convex_record_data1, convex_record_data1, proj4string=CRS("+proj=longlat +datum=WGS84"))
				convex_spatial1=spTransform(convex_spatial1, CRS("+proj=cea +lat_ts=30"))	
				convex_coordsUse1=convex_spatial1@coords
				colnames(convex_coordsUse1)<-c("X","Y")
				convex_convex1 = calcConvexHull(convex_coordsUse1)
				convex_area1=calcArea(convex_convex1)$area
				convex_area_km1<-convex_area1/10^6 # to get into square km.
			
				# save the result
				res_array_chull[k,j]<-convex_area_km1
				
				### then the grid cell method
				# raster should already be made...just need to make sure all components are in equal area projection
				# new projection
				albers<-"+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +el"

				# reconvert all the .shp files to equal area
				species_albers<-spTransform(c.species, CRS(albers))
				sites_albers<-spTransform(extinction_record, CRS(albers))
				stage_albers<-spTransform(c.period, CRS(albers))
				map_albers<-spTransform(US_map_rangeproj, CRS(albers))
				
				# and the raster
				proj4string(x.raster) <- CRS("+proj=longlat")
				raster_albers <- projectRaster(x.raster, crs=albers)
				
				# looking good - now extract various aspects...starting with baseline potential range
				# max.pot.occupancy<-extract(raster_albers, species_albers)
				# pot.grid.cells<-length(unique(max.pot.occupancy[[1]]))
				# print(pot.grid.cells)

				# actual occupancy of species
				actual.occupancy<-extract(raster_albers, sites_albers)
				actual.grid.cells<-length(unique(actual.occupancy))

				# expressed as a percentage:
				# grid.cell.perc<-actual.grid.cells/pot.grid.cells*100
				
				# save the results
				res_array_gcell[k,j]<-actual.grid.cells
				# res_array_gcell_perc[k,j]<-grid.cell.perc
						
				} # for all iterations

			} # only if the sedimentary record and species range overlap
			
			write.csv(res_array_latrange, lr.s.name, row.names=FALSE) # this is to save things, if the code breaks (which it always does)
			write.csv(res_array_gcd, gcd.s.name, row.names=FALSE)
			write.csv(res_array_chull, chull.s.name, row.names=FALSE)
			write.csv(res_array_gcell, gcell.s.name, row.names=FALSE)
			
		} # for each species	
		
	write.csv(res_array_latrange, lr.s.name, row.names=FALSE)
	write.csv(res_array_gcd, gcd.s.name, row.names=FALSE)
	write.csv(res_array_chull, chull.s.name, row.names=FALSE)
	write.csv(res_array_gcell, gcell.s.name, row.names=FALSE)
	# write.csv(res_array_gcell_perc, gcell.s.perc.name, row.names=FALSE)
		
	} # for each no. sites
	
	# stop the clock
	b<-toc()
	time_res[i,s]<-(b$toc/60/60)-(b$tic/60/60) # to get it into hours
	
} # for each period

time_res # just to give us a nice printout at the end of the run.
write.csv(time_res, "time_elapsed.csv", row.names=FALSE)

#####################################################################################################################
#####################################################################################################################

# small extra snippets to help with plotting, as a check to make sure it's doing what we want
us_map_clipped<-crop(US_map_rangeproj, extent(-135, -65, 22, 50))
us_map_clipped_albers<-spTransform(us_map_clipped, CRS(albers))

us_range<-gIntersection(us_map_clipped_albers, species_albers)

dev.new(height=8, width=10)
plot(us_map_clipped_albers, col="light grey", axes=T)
plot(us_range, col="light blue", add=T)
plot(stage_albers, add=T, col="yellow")
points(sites_albers, pch=19, col="red")

# trying to sort out the grid cell stuff
max.pot.occupancy<-extract(raster_albers, species_albers)
str.length<-length(max.pot.occupancy)
unique.cells.vector<-c()
for (shape in 1:str.length){
	this.shape<-as.vector(max.pot.occupancy[[shape]])
	unique.cells.vector<-c(unique.cells.vector, this.shape)
}
unique.cells.vector<-unique(unique.cells.vector)
pot.grid.cells<-length(unique.cells.vector)

# some other plotting stuff

t.species<-gIntersection(c.species, US_map_rangeproj)

# plot
dev.new(height=8, width=10)
plot(t.species, col="white", axes=T, main=sp.name, xlim=c(-103.001, -103.003), ylim=c(36.774, 36.776))
plot(US_map_rangeproj, col="light grey", add=T)
plot(t.species, col="light blue", add=T)
plot(species_cropped, add=T, col="yellow", border="yellow")
points(extinction_record, col="red", pch=19)

dev.new(height=8, width=10)
plot(t.species, col="white", axes=T, main=sp.name, xlim=c(-103.000, -103.005), ylim=c(36.772, 36.778))
plot(US_map_rangeproj, col="light grey", add=T)
plot(t.species, col="light blue", add=T)
plot(c.period, add=T, col="yellow", border="yellow")

require(devtools)
# install_github("eblondel/cleangeo")
require(cleangeo)

report<-clgeo_CollectionReport(Cretaceous) # shows that we have ring self-intersections...
summary <- clgeo_SummaryReport(report)
issues <- report[report$valid == FALSE,]
issues

#####################################################################################################################
#####################################################################################################################

### okay...try and set up a system whereby we can start the code at whatever species we want...

start.species <- "Sorex tenellus" # easy one to start with
remaining.sp.names <- sp.names[c(which(sp.names == start.species)) : length(sp.names)]
remainder_no <- "9" # remember to change this so you don't save over everything you've done

###

selected_NAtaxa<-list()
counter=1

for (t in 1:length(mam_unique)) {    
    ts<-mam_unique[t]
    tmp<-mam_data[which(mam_data$scientific == ts) , ]
    name<-as.character(tmp$scientific[1])
    # print(name)   
    if (name %in% remaining.sp.names == TRUE){    	
    	selected_NAtaxa[[counter]]<-tmp
    	counter=counter+1    	
    }  
    # print(counter)    
}   

###

remaining_species<-selected_NAtaxa
l.re.species<-length(remaining_species)

# get a vector species names
r.sp.names<-c()

for (smut in 1:length(selected_NAtaxa)){
	t.smut<-selected_NAtaxa[[smut]]
	nm<-t.smut$scientific[1]
	r.sp.names<-c(r.sp.names, as.character(nm))
}

# Analysis

for (i in 1:l.periods){
	
	per<-all_periods[[i]]
	c.period<-spTransform(per, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
	
	# identify the period and give ourselves a check
	p.name<-p.names[i]
	# print(p.name)
	
		for (s in 1:l.sites){
		
		# start the clock
		a<-tic()
		
		# establish no. sites	
		t.sites<-sites[s]
					
		### make some results arrays
		# make a results table for lat. range
		res_array_latrange<-array(NA, dim=c(iterations, l.re.species))
		colnames(res_array_latrange)<-c(as.character(r.sp.names))
		lr.s.name<-paste("remainder", remainder_no, p.name, "latrange", t.sites, "sites", "csv", sep=".")
	
		# make a results table for GCD
		res_array_gcd<-array(NA, dim=c(iterations, l.re.species))
		colnames(res_array_gcd)<-c(as.character(r.sp.names))
		gcd.s.name<-paste("remainder", remainder_no, p.name, "gcd", t.sites, "sites", "csv", sep=".")
		
		# make a results table for convex hull
		res_array_chull<-array(NA, dim=c(iterations, l.re.species))
		colnames(res_array_chull)<-c(as.character(r.sp.names))
		chull.s.name<-paste("remainder", remainder_no, p.name, "chull", t.sites, "sites", "csv", sep=".")
		
		# make a results table for grid cell
		res_array_gcell<-array(NA, dim=c(iterations, l.re.species))
		colnames(res_array_gcell)<-c(as.character(r.sp.names))
		gcell.s.name<-paste("remainder", remainder_no, p.name, "gcell", t.sites, "sites", grid_res, "deg", "csv", sep=".")
		
		# make a results table for grid cell percentage - not worrying about this for now
		# res_array_gcell_perc <-array(NA, dim=c(iterations, l.species))
		# colnames(res_array_gcell_perc)<-c(as.character(sp.names))
		# gcell.s.perc.name<-paste(p.name, "gcell", t.sites, "sites", grid_res, "percentage", "deg", "csv", sep=".")
	
			for (j in 1: l.re.species){
		
			c.species<-remaining_species[[j]]
			projection(c.species)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"	# get the range in the right projection
				
			# cut the polygons
			species_cropped <-gIntersection(c.period, c.species)
		
			sp.name<-r.sp.names[j]
			# print(sp.name)
			
			# print a check string
			check<-paste(p.name, t.sites, sp.name, sep="    ")
			print(check)
		
			if (is.null(species_cropped)) { 			
				
				res_array_latrange[,j] <- c(NA)
				res_array_gcd[,j] <- c(NA)
				res_array_chull[,j] <- c(NA)
				res_array_gcell[,j] <- c(NA)
			
			} else {		
	
				for (k in 1:iterations){
			
				# print(k)
				
				# simulate a fossil record
				extinction_record<-spsample(species_cropped, sites, type="random", iter=1000000)
			
				### start with lat. range
				latlong_record<-cbind(extinction_record$x, extinction_record$y)
				max.lat<-max(latlong_record[,2])
				min.lat<-min(latlong_record[,2])
				sim.lat.range<-abs(max.lat - min.lat)
			
				# save the result
				res_array_latrange[k,j]<-sim.lat.range
				
				### then GCD
				# set up an empty vector and organize the coordinates
				gcdists<-vector()
				gc_record<-cbind(extinction_record$x, extinction_record$y)
		
				# then create a nested loop to perform a great circle cal. (spDistsN1) on each point, and tape the results to the vector
				for (row in 1:nrow(gc_record)){
					now.row<-gc_record[row,]
					gcs<-spDistsN1(gc_record, now.row, longlat = TRUE)
					gcdists<-c(gcdists, gcs)						
				}

				# now need to get the maximum of all these possible combinations, and work out how close it is to the actual maximum great circle distance of the range
				max.gcdist<-max(gcdists)
				res_array_gcd[k,j]<-max.gcdist
				
				### then convex hull
				convex_record_alt1<-cbind(extinction_record$x, extinction_record$y)
				convex_record_data1<-as.data.frame(convex_record_alt1)
				convex_spatial1<-SpatialPointsDataFrame(convex_record_data1, convex_record_data1, proj4string=CRS("+proj=longlat +datum=WGS84"))
				convex_spatial1=spTransform(convex_spatial1, CRS("+proj=cea +lat_ts=30"))	
				convex_coordsUse1=convex_spatial1@coords
				colnames(convex_coordsUse1)<-c("X","Y")
				convex_convex1 = calcConvexHull(convex_coordsUse1)
				convex_area1=calcArea(convex_convex1)$area
				convex_area_km1<-convex_area1/10^6 # to get into square km.
			
				# save the result
				res_array_chull[k,j]<-convex_area_km1
				
				### then the grid cell method
				# raster should already be made...just need to make sure all components are in equal area projection
				# new projection
				albers<-"+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +el"

				# reconvert all the .shp files to equal area
				species_albers<-spTransform(c.species, CRS(albers))
				sites_albers<-spTransform(extinction_record, CRS(albers))
				stage_albers<-spTransform(c.period, CRS(albers))
				map_albers<-spTransform(US_map_rangeproj, CRS(albers))
				
				# and the raster
				proj4string(x.raster) <- CRS("+proj=longlat")
				raster_albers <- projectRaster(x.raster, crs=albers)
				
				# looking good - now extract various aspects...starting with baseline potential range - not going to do this for now (easier ways of doing it outside loop)
				# max.pot.occupancy<-extract(raster_albers, species_albers)
				# str.length<-length(max.pot.occupancy)
				# unique.cells.vector<-c() # we have to do this because the 'extract' function is giving us a list of grid cells occupied by each polygon
				# for (shape in 1:str.length){
				#	this.shape<-as.vector(max.pot.occupancy[[shape]])
				#	unique.cells.vector<-c(unique.cells.vector, this.shape)
				#	}
				# unique.cells.vector<-unique(unique.cells.vector)
				# pot.grid.cells<-length(unique.cells.vector)

				# actual occupancy of species
				actual.occupancy<-extract(raster_albers, sites_albers)
				actual.grid.cells<-length(unique(actual.occupancy))

				# expressed as a percentage:
				# grid.cell.perc<-actual.grid.cells/pot.grid.cells*100
				
				# save the results
				res_array_gcell[k,j]<-actual.grid.cells
				# res_array_gcell_perc[k,j]<-grid.cell.perc
						
				} # for all iterations

			} # only if the sedimentary record and species range overlap
			
			write.csv(res_array_latrange, lr.s.name, row.names=FALSE) # this is to save things, if the code breaks (which it always does)
			write.csv(res_array_gcd, gcd.s.name, row.names=FALSE)
			write.csv(res_array_chull, chull.s.name, row.names=FALSE)
			write.csv(res_array_gcell, gcell.s.name, row.names=FALSE)	
			
		} # for each species	
		
	write.csv(res_array_latrange, lr.s.name, row.names=FALSE)
	write.csv(res_array_gcd, gcd.s.name, row.names=FALSE)
	write.csv(res_array_chull, chull.s.name, row.names=FALSE)
	write.csv(res_array_gcell, gcell.s.name, row.names=FALSE)
	# write.csv(res_array_gcell_perc, gcell.s.perc.name, row.names=FALSE)
		
	} # for each no. sites
	
	# stop the clock
	b<-toc()
	time_res[i,s]<-(b$toc/60/60)-(b$tic/60/60) # to get it into hours
	
} # for each period

time_res # just to give us a nice printout at the end of the run.
time.title<-paste(p.name, t.sites, "time_elapsed", "csv", sep=".")
write.csv(time_res, time.title, row.names=FALSE)

#####################################################################################################################
#####################################################################################################################

# works great

### screwing about with plotting and aggregation
species<-all_species[[20]]

# experiment with site aggregation
extinction_record_1<-spsample(species, 20, type="random", iter=1000000)
extinction_record_2<-spsample(species, 20, type="clustered", nclusters=3, iter=1000000)
extinction_record_1
extinction_record_2

dev.new(width=6, height=8)
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))

plot(US_map_rangeproj, col="light grey", xlim=c(-130, -65), ylim=c(20, 50), axes=T)
plot(species, col="white", add=T)
points(extinction_record_1, pch=21, col="black", bg="red")

plot(US_map_rangeproj, col="light grey", xlim=c(-130, -65), ylim=c(20, 50), axes=T)
plot(species, col="white", add=T)
points(extinction_record_2, pch=21, col="black", bg="red")

### look in a loop
dev.new(width=10, height=8)
par(mfrow=c(3,3))
par(mar=c(2,2,2,2))

for (t in 1:9){
	
	extinction_record_2<-spsample(species, 20, type="clustered", nclusters=4, iter=1000000)
	print(extinction_record_2)
	plot(US_map_rangeproj, col="light grey", xlim=c(-130, -65), ylim=c(20, 50), axes=T)
	plot(species, col="white", add=T)
	points(extinction_record_2, pch=21, col="black", bg="red")
	
}


### not working brilliantly, try package 'spatstat'
projection(species)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projection(Cretaceous)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
species_cropped<-gIntersection(Cretaceous, species)
range_cropped<-gIntersection(US_map_rangeproj, species)

dev.new(width=10, height=8)
par(mfrow=c(3,3))
par(mar=c(2,2,2,2))

for (t in 1:9){
	
	extinction_record_1<-spsample(species_cropped, 20, type="random", iter=1000000)
	extinction_record_2<-spsample(species_cropped, 20, type="clustered", nclusters=4, iter=1000000)
	print(extinction_record_2)
	plot(US_map_rangeproj, col="light grey", xlim=c(-130, -65), ylim=c(20, 50), axes=T)
	plot(range_cropped, col="white", add=T)
	points(extinction_record_2, pch=21, col="black", bg="red")
	
}

extinction_record_1<-spsample(species_cropped, 20, type="random", iter=1000000)
extinction_record_2<-spsample(species_cropped, 20, type="clustered", nclusters=4, iter=1000000)

# okay, now work in cutting out Cretaceous rock
dev.new(width=10, height=8)
par(mfrow=c(1,1))
par(mar=c(2,2,2,2))

plot(US_map_rangeproj, col="light grey", xlim=c(-130, -65), ylim=c(20, 50), axes=T)
plot(range_cropped, col="white", add=T)
plot(species_cropped, col="light blue", border="light blue", add=T)
points(extinction_record_2, pch=21, col="black", bg="red")



