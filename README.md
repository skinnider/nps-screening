## NPS screening

This repository contains code used for retrospective analysis of full-scan mass spectrometry data to prioritize emerging NPS for assay development. 

The workflow is as follows:
- `R/MS/outer-write-ms1-matches-xcms.R`: searches a series of mzML mass spectrometry files one at a time for fragment ions from a suspect database, retaining spectra that match at least two fragment ions from a database entry. Designed to be run on a computing platform that supports array jobs (e.g., slurm, Torque, PBSPro) using functions in the `R/functions` directory.
- `R/MS/outer-consolidate-ms1-matches-xcms.R`: consolidates the results from `write-ms1-matches-xcms` into a single spreadsheet, in RDS format. 
- `R/figures/plot-fig*.R`: code used to generate the visualizations shown in the manuscript. 

Key data files are as follows:

- `NPS DATABASE-NOV2022.csv`: this is the suspect database containing published fragment ions for 83 NPS analyzed in the study.
- `data/databases/DRUGS-200323.csv`: MS/MS spectra collected at the PTC from newly-acquired reference standards. 
- `data/xcms/snthresh=10-noise=100-ppm=25-peakwidth_min=5-peakwidth_max=20.rds`: this is the RDS file output by `consolidate-ms1-matches-xcms`. It contains patient information so is not provided in this repository, but any code referencing this file should also work with a newly-generated RDS file. 

The code was run in a conda environment that was created using the code in `sh/create-env.sh`. 
