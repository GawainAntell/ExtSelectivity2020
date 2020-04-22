# Extinction selectivity project overview
This repository contains all data and R code necessary to evaluate extinction selectivity on species' geographic range size. The outputs of the scripts are provided in the /Results folder. Please cite the study as:

Darroch, S.D., Casey, M., Antell, G.S., Sweeney, A., and Saupe, E.E. (in press). High preservation potential of paleogeographic range size distributions in deep time. American Naturalist.

Simulate fossil species distributions from by-sediment site placement:
- convert_mam_polygons_to_brick_GSA.R
- simulate_by_sed_GSA.R

Simulate fossil species distributions from by-species site placement:
- contiguous_US_code_clean_SD.R
- reformat_SD_data.R

Note that the simulation files require raw data inputs that are not supplied in this repository, including raster files exported from ArcGIS. These files will be archived at Zenodo in zipped format. The simulation files listed above manipulate the input data and export 5 intermediate data files that are then called by the analysis scripts.

Analyze accuracy of paleo-range size reconstruction:
- tau_analysis_GSA.R

Analyze accuracy and precision of extinction selectivity:
- beta_analysis_GSA.R


The data necessary to source the 2 analysis scripts are provided in the /Data folder:
- DS1_stage_data.csv: age and coverage information for the 13 Carboniferous stages
- DS2_IUCN_range_data_vector.csv: true geographic range sizes measured directly from polygons
- DS3_IUCN_range_data_0.5degree_raster.csv: true geographic range sizes measured after rasterization
- DS4_simulation_range_data_by_species.csv: simulated fossil distributions, by-species site placement
- DS5_simulation_range_data_0.5_degree_by_sed.csv: simulated fossil distributions, by-sediment site placement
