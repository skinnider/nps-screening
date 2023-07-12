# Find chromatographic peaks with xcms and extract ones with a precursor mass
# and at least one fragment matching a compound in the NPS database.
setwd("~/git/nps-screening")
options(stringsAsFactors = FALSE)
library(argparse)

# dynamically parse arguments
parser = ArgumentParser(prog = 'inner-write-ms1-matches-xcms.R')
grid = read.delim("sh/grids/write-ms1-matches-xcms.txt")
for (param_name in colnames(grid))
  parser$add_argument(paste0('--', param_name),
                      type = typeof(grid[[param_name]]))
args = parser$parse_args()
print(args)

library(tidyverse)
library(magrittr)
library(xcms)

# create metadata frame
meta = data.frame(file = args$mzml_file)

# read data
dat = readMSData(files = args$mzml_file, 
                 pdata = new("NAnnotatedDataFrame", meta),
                 mode = "onDisk")

# run peak detection
## can probably run with a grid of snthresh, peakwidth, ppm, noise
cwp = CentWaveParam(snthresh = args$snthresh,
                    noise = args$noise, 
                    ppm = args$ppm,
                    peakwidth = c(args$peakwidth_min, args$peakwidth_max))
dat = findChromPeaks(dat, param = cwp)

# extract MS/MS spectra
spectra = chromPeakSpectra(dat, msLevel = 2L)

# save all MS/MS spectra detected in this file
spectra_df = data.frame(spectrum_name = names(spectra)) %>% 
  separate(spectrum_name, into = c('chrom_peak', 'file', 'spectrum'),
           sep = "\\.", remove = FALSE) %>% 
  dplyr::select(-file)

# now, read databases
db1 = read.csv("data/databases/NPS DATABASE-NOV2022.csv", skip = 5)
databases = list(NPS = db1) %>% 
  # filter to MS1/MS2 rows only
  map(~ {
    db = .x
    stop_at = which(db$Compound.Name == "") %>% head(1)
    db %<>% extract(seq_len(stop_at), )
    keep = map_lgl(db, ~ n_distinct(.x) > 1)
    db %<>% extract(, keep)
  })

# function to calculate ppm boundary
calc_ppm_range = function(theor_mass, err_ppm = 10) {
  c(
    (-err_ppm / 1e6 * theor_mass) + theor_mass,
    (err_ppm / 1e6 * theor_mass) + theor_mass
  )
}

# iterate through compounds
compounds = unique(databases$NPS$Compound.Name) %>% na.omit()
## remove one compound already in the in-house database
compounds %<>% setdiff(c('', 'Norfluorodiazepam'))
results = map(seq_along(compounds), ~ {
  compound = compounds[[.x]]
  message("[", .x, "/", length(compounds), "] ", compound, " ...")
  
  # get parent compound info to extract candidate spectra
  compound = filter(databases$NPS, Compound.Name == compound)
  parent = filter(compound, Workflow == 'TargetPeak')
  fragments = filter(compound, Workflow == 'Fragment')
  mz_range = calc_ppm_range(parent$m.z, err_ppm = 10)

  # find spectra that match parent properties
  ms1_match = map_lgl(seq_along(spectra), ~ {
    spectrum = spectra[[.x]]
    between(precursorMz(spectrum), mz_range[1], mz_range[2]) & 
      precursorIntensity(spectrum) >= parent$Height.Threshold
  })
  
  # for the spectra that match, filter to those that contain a matching fragment
  ms2_match = map_lgl(seq_along(spectra), ~ {
    if (!ms1_match[.x]) return(FALSE)
    spectrum_df = spectra[[.x]] %>% as.data.frame()
    match = tidyr::crossing(spectrum_df, target_mz = fragments$Product.m.z) %>% 
      mutate(match = map2_lgl(mz, target_mz, ~ {
        range = calc_ppm_range(.y, err_ppm = 20)
        between(.x, range[1], range[2])
      }))
    any(match$match)
  })
  
  # abort if no matches at all
  if (!any(ms2_match)) return(list())
  
  # extract just those spectra
  spectrum_names = names(spectra)[ms2_match]
  spectrum_ms1 = map(which(ms2_match), ~ {
    spectrum = spectra[[.x]]
    data.frame(mz = precursorMz(spectrum),
               rt = rtime(spectrum),
               intens = precursorIntensity(spectrum))
  }) %>% 
    setNames(spectrum_names) %>% 
    bind_rows(.id = 'spectrum')
  spectrum_ms2 = map(which(ms2_match), ~ as.data.frame(spectra[[.x]])) %>% 
    setNames(spectrum_names) %>% 
    bind_rows(.id = 'spectrum')
  output = list(ms1 = spectrum_ms1, ms2 = spectrum_ms2) %>% 
    map(~ cbind(parent, .x))
}) %>% 
  setNames(compounds)

# add chromPeakSpectra results
output = list(chromPeakSpectra = spectra_df,
              ms1_matches = results)

# save
output_dir = dirname(args$output_file)
if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
saveRDS(output, args$output_file)
