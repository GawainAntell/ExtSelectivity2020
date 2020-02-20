# ExtSelectivity
This repository contains all data and R code necessary to evaluate extinction selectivity on species' geographic range size, as implemented in Darroch et al. 2020. The outputs of the scripts are provided in the /Results folder.

Simulate fossil species distributions from by-sediment site placement:
convert_mam_polygons_to_brick_GSA.R
simulate_by_sed_GSA.R
Simulate fossil species distributions from by-species site placement:
contiguous_US_code_clean_SD.R
reformat_SD_data.R
Analyze accuracy of paleo-range size reconstruction:
tau_analysis_GSA.R
Analyze accuracy and precision of extinction selectivity:
beta_analysis_GSA.R

The data necessary to source the analysis scripts are provided in the /Data folder:
DS1_stage_data.csv: age and coverage information for the 13 Carboniferous stages
DS2_IUCN_range_data_vector.csv: true geographic range sizes measured directly from polygons
DS3_IUCN_range_data_0.5degree_raster.csv: true geographic range sizes measured after rasterization
DS4_simulation_range_data_by_species.csv: simulated fossil distributions, by-species site placement
DS5_simulation_range_data_0.5_degree_by_sed.csv: simulated fossil distributions, by-sediment site placement
