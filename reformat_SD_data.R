library(stringr)
# This script compiles the dozens of files output from running
# contiguous_US_code_clean_SD.R, that is, the by-species simulation data.
# The compiled csv files are named DS1 and DS2, and these are
# saved in the Data folder for use in the tau and beta analysis scripts.
# To run this code, first unzip the 'SD_data' folder.

foldrExpt <- 'Data/data_by_species/'
foldrMain <- 'Data/SD_data/'  
foldrStage <- paste0(foldrMain,'For_multivariate/')

####################################
# extract true range of all species
iucn <- read.csv(paste0(foldrMain,'baseline.res.csv'), stringsAsFactors=FALSE)
rangeCols <- which(colnames(iucn) %in% c("area.km2.","GCD","gridcells","lat.range"))
newNames <- c('chull','latrange','gcd','gcell')
colnames(iucn)[rangeCols] <- newNames
colnames(iucn) <- tolower(colnames(iucn))
noRange <- is.na(iucn$chull)
iucn <- iucn[!noRange,]
iucn$species <- gsub(' ', '_', iucn$species)
spOrdr <- order(iucn$species)
iucn <- iucn[spOrdr,]
iucn2print <- iucn[,c('species',newNames)]
write.csv(iucn2print, paste0(foldrExpt,'DS2_IUCN_range_data_vector.csv'), row.names = FALSE)

####################################
# collate fossil range area for each given metric and number of sites

# each file contains 100 iterations of range size for a metric, stage & n sites
foldrFslRange <- paste0(foldrMain, 'results_csv_files_carboniferous/')
toRead <- list.files(foldrFslRange)

# remove metadata files
toRead <- toRead[-grep('elapsed', toRead)]

exptFsl <- stageSppCount <- data.frame()
for (nm in toRead){
  fsl <- read.csv(paste0(foldrFslRange, nm)) 
  fsl <- data.frame( t(fsl) )
  
  iterColNms <- paste0('rep',1:100)
  colnames(fsl) <- iterColNms
  
  # format fossil data to match IUCN structure
  fslSpOrdr <- order(rownames(fsl))
  fsl <- fsl[fslSpOrdr,]
  fsl$species <- gsub('\\.','_', row.names(fsl))
  
  # note there are more species in the fossil dataset than have true range size
  # there are also unpreserved species in the fossil dataset (range = NA)
  sameSp <- fsl$species %in% iucn$species
  fsl <- fsl[sameSp,]
  prsrvd <- !is.na(fsl$rep1)
  fsl <- fsl[prsrvd,]
  
  # include simulation attribute data
  params2add <- unlist(str_split(nm, '\\.'))[1:3]
  atrM <- matrix(rep(params2add, nrow(fsl)), ncol=3, byrow = TRUE)
  atrCols <- c('stage','metric','sites')
  fsl[,atrCols] <- atrM
  stageAbrev <- substr(params2add[1], 1, 2)
  s <- sprintf('%03d', as.numeric(params2add[3]))
  m <- params2add[2]
  fsl$sim_id <- paste0(stageAbrev, s, m)
  fsl <- fsl[,c('species',atrCols, iterColNms,'sim_id')]
  
  exptFsl <- rbind(exptFsl, fsl)
  
  # save number of preserved species; original file was inaccurate
  nSpp <- sum(prsrvd)
  stageSppCount <- rbind(stageSppCount, cbind(params2add[1],nSpp))
}

exptFsl$stage <- gsub('Guadelupian','Guadalupian',exptFsl$stage)
write.csv(exptFsl, paste0(foldrExpt,'DS4_simulation_range_data_by_species.csv'), row.names=FALSE)

colnames(stageSppCount) <- c('stage','nSpp')
keep <- ! duplicated(stageSppCount$stage)
stageSppCount <- stageSppCount[keep,]
stageSppCount$stage <- as.character(stageSppCount$stage)
typo <- stageSppCount$stage=='Guadelupian'
stageSppCount$stage[typo] <- 'Guadalupian'
countNm <- paste0(foldrStage,'no.sp.stage_new.csv')
write.csv(stageSppCount, countNm, row.names = FALSE)

####################################
# collate metadata on stages

stageA <- read.csv(paste0(foldrStage, 'carb.areas.csv'), stringsAsFactors=FALSE)
stageMst <- read.csv(paste0(foldrStage, 'min.span.trees.stage.csv'), stringsAsFactors=FALSE)
stageAge <- read.csv(paste0(foldrStage, 'NAm_stage_chronology_GSA.csv'), stringsAsFactors=FALSE)
stageSppCount <- read.csv(paste0(foldrStage, 'no.sp.stage_new.csv'), stringsAsFactors=FALSE)
colnames(stageA)[-1] <- colnames(stageMst)[-1] <- colnames(stageSppCount) <- c('stage','value')
nmOrdr <- match(stageAge$stage, stageSppCount$stage)
stageSppCount <- stageSppCount[nmOrdr,]

stageDat <- cbind(stageAge[,c('stage','duration','age_mid')], 
                  stageA$value, stageMst$value, stageSppCount$value)
colnames(stageDat) <- c('stage','duration','midpoint','area','dispersion','speciesN')
write.csv(stageDat, paste0(foldrExpt,'DS1_stage_data.csv'), row.names=FALSE)
