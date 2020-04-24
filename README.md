# Extinction selectivity project overview

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3765445.svg)](https://doi.org/10.5281/zenodo.3765445)

This repository contains all data and R code necessary to evaluate extinction selectivity on species' geographic range size. The outputs of the scripts are provided in the /Results folder. Please cite the study as:

Darroch, S.A.F., Casey, M.M., Antell, G.S., Sweeney, A., and Saupe, E.E. (2020). High preservation potential of paleogeographic range size distributions in deep time. American Naturalist.

There are 4 R files that manipulate the raw input data and export the 5 intermediate data files that are then called by the analysis scripts. Before running the simulation scripts, download the zipped /spatial_data and /SD_data folders at https://doi.org/10.5281/zenodo.3765445 and move them to /Data. Note that due to the probabilistic nature of some simulation steps (e.g. the selection of species to survive vs. go extinct), the exact output of each simulation will vary slightly from one run to another.

Simulate fossil species distributions from by-sediment site placement:
- convert_mam_polygons_to_brick_GSA.R
- simulate_by_sed_GSA.R

Simulate fossil species distributions from by-species site placement:
- contiguous_US_code_clean_SD.R
- reformat_SD_data.R

There are 2 R files that produce the final figures and tables from the intermediate data files (DS1 - DS5, described below). Since these intermediate files are provided, it is not necessary to have run the simulation code files listed above.

Analyze accuracy of paleo-range size reconstruction:
- tau_analysis_GSA.R

Analyze accuracy and precision of extinction selectivity:
- beta_analysis_GSA.R

The following intermediate data files, which are output by the simulation code, and which are input to the analysis scripts, are provided in the /Data folder:
- DS1_stage_data.csv: age and coverage information for the 13 Carboniferous stages
- DS2_IUCN_range_data_vector.csv: true geographic range sizes measured directly from polygons
- DS3_IUCN_range_data_0.5degree_raster.csv: true geographic range sizes measured after rasterization
- DS4_simulation_range_data_by_species.csv: simulated fossil distributions, by-species site placement
- DS5_simulation_range_data_0.5_degree_by_sed.csv: simulated fossil distributions, by-sediment site placement

Please report any issues with code or data files either in the GitHub repo (https://github.com/GwenAntell/ExtSelectivity2020) or to the corresponding author (simon.darroch@gmail.com).
