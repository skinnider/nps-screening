setwd("~/git/nps-screening")
options(stringsAsFactors = FALSE)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-consolidate-ms1-matches-xcms.R')
grid = read.delim("sh/grids/consolidate-ms1-matches-xcms.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(xcms)

source("R/functions/detect_system.R")

# define input directory
input_dir = file.path(base_dir, "UDS", "converted") %>%
  gsub("scratch", "project", .)

# list mzML files
mzml_files = list.files(input_dir, full.names = TRUE, pattern = "mzML") %>% 
  # drop calibration, mix, QC files
  extract(!grepl('^Cal|^MeOH|^QC|^Mix', basename(.))) %>% 
  # keep only E/R 
  extract(grepl('^E|^R', basename(.)))

# get peak-picked files
mzml_filenames = gsub("\\.mzML$", "", basename(mzml_files))
rds_files = file.path(args$output_dir, paste0(mzml_filenames, '.rds'))
stopifnot(all(file.exists(rds_files)))

# read all files
dats = map(rds_files, readRDS) %>% setNames(mzml_filenames)
# extract spectra
chromPeakSpectra = map(dats, 'chromPeakSpectra') %>% 
  bind_rows(.id = 'file')
# extract matches
matches_ms1 = map(dats, 'ms1_matches') %>% 
  map(~ extract(., lengths(.) > 0) %>% map('ms1') %>% bind_rows()) %>% 
  bind_rows(.id = 'file')
matches_ms2 = map(dats, 'ms1_matches') %>% 
  map(~ extract(., lengths(.) > 0) %>% map('ms2') %>% bind_rows()) %>% 
  bind_rows(.id = 'file')

# save all to RDS
output = list(chromPeakSpectra = chromPeakSpectra,
              matches = list(ms1 = matches_ms1, ms2 = matches_ms2))

# save
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(output, args$output_file)

# remove all other RDS files
if (file.exists(args$output_file))  {
  dat = readRDS(args$output_file)
  if (length(dat) == 2 & 
      all(names(dat) %in% c('chromPeakSpectra', 'matches')))
    file.remove(rds_files)
}
