setwd("~/git/nps-screening")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(Spectra)
source("R/theme.R")

# read standards
stds = read.csv("data/databases/DRUGS-200323.csv", skip = 2)
# convert to spectra
std_df = pmap_dfr(stds, function(...) {
  current = tibble(...)
  parent = current$ExtractedMass
  fragments = dplyr::select(current, starts_with('Fragment')) %>% 
    map_dbl(identity)
  intens = dplyr::select(current, starts_with('Target.Ratio')) %>% 
    map_dbl(identity)
  df = data.frame(mz = fragments, intensity = intens) %>% 
    mutate(compound = current$CompoundName, parent = parent) %>% 
    dplyr::select(compound, parent, mz, intensity) %>% 
    arrange(mz) %>% 
    # normalize fragment intensities
    mutate(intensity = intensity / max(intensity)) %>% 
    # remove rownames
    remove_rownames()
})
std_spd = std_df %>% 
  distinct(compound, parent) %>% 
  mutate(msLevel = 2L, rtime = 1) %>% 
  arrange(compound)
std_spd$mz = std_df %>% split(.$compound) %>% map('mz')
std_spd$intensity = std_df %>% split(.$compound) %>% map('intensity')
stopifnot(all(std_spd$spectrum == names(std_spd$mz))) ## TRUE
stopifnot(all(std_spd$spectrum == names(std_spd$intensity))) ## TRUE
std_sps = Spectra(std_spd)
std_sps$name = std_spd$compound

# read all matches
dat = readRDS("data/xcms/snthresh=10-noise=100-ppm=25-peakwidth_min=5-peakwidth_max=20.rds")
ms1 = dat$matches$ms1
ms2 = dat$matches$ms2
# map to dates (ignore QC, etc.)
dates = readRDS("data/MS/dates.rds") %>% 
  dplyr::rename(file = patient) %>% 
  dplyr::select(file, date) %>% 
  arrange(date) %>% 
  group_by(file) %>% 
  dplyr::slice(1) %>% 
  ungroup()
ms1 %<>% left_join(dates) %>% drop_na(date)
ms2 %<>% left_join(dates) %>% drop_na(date)

# filter to compounds in the validation set
ms1 %<>% filter(Compound.Name %in% std_df$compound)
ms2 %<>% filter(Compound.Name %in% std_df$compound)

# now, read database
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

# calculate level of evidence for each peak
fragments = databases$NPS %>% 
  filter(MS.Order == 'ms2') %>% 
  distinct(Compound.Name, Product.m.z)
evidence_per_fragment = ms2 %>% 
  dplyr::select(-Product.m.z) %>% 
  left_join(fragments, by = 'Compound.Name') %>% 
  ## low-res filter
  filter(abs(mz - Product.m.z) < 1) %>% 
  group_by(file, Compound.Name, spectrum, mz) %>% 
  summarise(match = map2_lgl(mz, Product.m.z, ~ {
    range = calc_ppm_range(.y, err_ppm = 20)
    between(.x, range[1], range[2])
  }) %>% sum()) %>% 
  ungroup()
evidence_per_parent = evidence_per_fragment %>% 
  group_by(file, Compound.Name, spectrum) %>% 
  summarise(matches = sum(match)) %>% 
  ungroup() 
# also calculate proportion of reference product ions detected
total_fragments = dplyr::count(fragments, Compound.Name,
                               name = 'total_fragments')
evidence_per_parent %<>% left_join(total_fragments, by = 'Compound.Name')

# merge into ms1 data frame
ms1 %<>% left_join(evidence_per_parent, by = c('file', 'Compound.Name', 
                                               'spectrum'))

# filter to spectra with at least two matches
ms1_subset = filter(ms1, matches >= 2)

# compute dot-products for each
msms = ms2 %>% 
  dplyr::select(-(Compound.Name:Integration.Strategy)) %>% 
  inner_join(ms1_subset, by = c('file', 'spectrum', 'date'))
dot_product = map_dfr(unique(std_df$compound), ~ {
  compound = .x
  message(compound)
  msms0 = filter(msms, Compound.Name == compound) %>% 
    # rename fragment ion m/z
    dplyr::rename(mz = mz.x) %>% 
    # normalize fragment intensities
    group_by(file, spectrum, date) %>% 
    mutate(intens = i / max(i)) %>% 
    ungroup() %>% 
    # paste file into spectrum name
    mutate(spectrum_ = paste(file, spectrum))
  
  ## create spectra
  spd = distinct(msms0, spectrum_, rt) %>% 
    mutate(msLevel = 2L) %>% 
    dplyr::rename(rtime = rt) %>% 
    arrange(spectrum_)
  spd$mz = msms0 %>% split(.$spectrum_) %>% map('mz')
  spd$intensity = msms0 %>% split(.$spectrum_) %>% map('intens')
  stopifnot(all(spd$spectrum_ == names(spd$mz))) ## TRUE
  stopifnot(all(spd$spectrum_ == names(spd$intensity))) ## TRUE
  sps = Spectra(spd)
  sps$name = spd$spectrum_
  ## skip if <2 hits
  if (length(sps) < 2) return(data.frame())
  
  # calculate dot-products
  dp = compareSpectra(std_sps, sps, fun = 'dotproduct', ppm = 20) %>% 
    set_rownames(std_sps$name) %>% 
    set_colnames(sps$name)
  # convert to data frame
  df = reshape2::melt(dp, varnames = c('reference_spectrum', 'uds_spectrum'),
                      value.name = 'dot_product') %>% 
    mutate(reference_spectrum = std_sps$name[reference_spectrum] %>% as.character(),
           uds_spectrum = sps$name[uds_spectrum] %>% as.character()) %>% 
    mutate(compound = compound)
  return(df)
})

# keep only best dot-product per spectrum
best = dot_product %>%
  mutate(patient = gsub(" CP.*$", "", uds_spectrum),
         uds_spectrum = gsub("^.*CP", "", uds_spectrum)) %>% 
  group_by(patient, reference_spectrum) %>% 
  arrange(desc(dot_product)) %>% 
  dplyr::slice(1) %>% 
  ungroup()

###############################################################################-
## Boxplot: dot-products to reference spectra
###############################################################################-

pal = c(pals::stepped()[seq(1, 20, 4)],
        pals::stepped3()[seq(1, 20, 4)])
p1 = dot_product %>% 
  filter(reference_spectrum == compound) %>% 
  ggplot(aes(x = reference_spectrum, y = dot_product, fill = reference_spectrum,
             color = reference_spectrum)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.6, size = 0.35) +
  geom_jitter(height = 0, width = 0.15, shape = 21, size = 0.4, stroke = 0.2,
              fill = NA) +
  scale_x_discrete('Reference spectrum') +
  scale_y_continuous('Dot product', limits = c(0, 1)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1.5,
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))
p1
ggsave("fig/final/figure5a.pdf", p1, width = 5.5, height = 5.5, 
       units = "cm", useDingbats = FALSE)

# print means
library(broom)
dot_product %>% 
  filter(reference_spectrum == compound) %>% 
  group_by(compound) %>% 
  do(tidy(summary(.$dot_product)))

###############################################################################-
## Boxplot: dot-products to reference spectra
###############################################################################-

# pick spectra
manual_spectra = c(
  'Eutylone' = 'E3510208159',
  'Fluorofentanyl' = 'R3530154883',
  'Bromazolam' = 'R3550062527',
  'Furanyl UF-17' = 'R3370209295',
  'Deschloroetizolam' = 'E3490168435',
  ## unusually bad/good spectra
  'Eutylone' = 'R3310075610',
  'Deschloroetizolam' = 'R3540210870'
) %>%
  data.frame(Compound.Name = names(.), file = .)
sample = inner_join(ms1, manual_spectra) %>% 
  # keep best match per file
  filter(spectrum %in% best2$spectrum) 
distinct(sample, file, Compound.Name, spectrum) 

# plot MS/MS for the sample
msms = ms2 %>% 
  dplyr::select(-(Compound.Name:Integration.Strategy)) %>% 
  inner_join(sample, by = c('file', 'spectrum', 'date')) 
plots = pmap(manual_spectra, ~ {
  current = tibble(...)
  msms0 = inner_join(msms, current) %>% 
    # rename fragment ion m/z
    dplyr::rename(mz = mz.x) %>% 
    group_by(file, spectrum, date) %>% 
    mutate(intens = i / max(i)) %>%
    ungroup()
  
  # also get dot-product
  dp = filter(dot_product, reference_spectrum == current$Compound.Name, 
              compound == current$Compound.Name, 
              gsub(" CP.*$", "", uds_spectrum) == current$file)
  label = paste0('dot-product = ',
                 formatC(dp$dot_product, digits = 2, format = 'f'))
  
  # get reference
  ref = std_df %>% 
    filter(compound == current$Compound.Name) %>% 
    mutate(intens = -1 * intensity)
  
  # create mirror data
  mirror_df = ref %>% 
    tidyr::crossing(distinct(ms2_sample0, file, spectrum, date)) %>% 
    mutate(label = 'Reference') %>% 
    bind_rows(ms2_sample0 %>% mutate(label = 'Patient'))
  # label facets
  labels = distinct(mirror_df, file, facet) %>% 
    drop_na(facet) %$%
    setNames(facet, file)
  
  # flag matched fragments
  ref$match = map_lgl(ref$mz, ~ {
    range = calc_ppm_range(.x, err_ppm = 20)
    map_lgl(msms0$mz, ~ between(.x, range[1], range[2])) %>% any()
  })
  msms0$match = map_lgl(msms0$mz, ~ {
    range = calc_ppm_range(.x, err_ppm = 20)
    map_lgl(ref$mz, ~ between(.x, range[1], range[2])) %>% any()
  })
  
  # create mirror data
  mirror_df = ref %>% 
    tidyr::crossing(distinct(msms0, file, spectrum, date)) %>% 
    mutate(label = 'Database') %>% 
    bind_rows(msms0 %>% mutate(label = 'Patient'))
  
  pal = pals::cols25()
  pal = c('FALSE' = 'black', 'TRUE' = pals::cols25()[2])
  bg = distinct(mirror_df, file, spectrum, date)
  p1 = mirror_df %>%
    ggplot() +
    ggtitle(current$Compound.Name, label) +
    geom_blank() + 
    geom_rect(data = bg,
              aes(xmin = -Inf, xmax = Inf, ymax = 0, ymin = -Inf),
              alpha = 0.1, fill = 'grey70', color = NA) +
    geom_label(data = bg %>% mutate(label = 'Reference'),
               aes(x = Inf, y = -Inf, label = label), hjust = 1, vjust = 0,
               size = 1.75, label.size = NA, label.padding = unit(0.5, 'lines'),
               fill = NA) +
    geom_label(data = bg %>% mutate(label = 'UDS'),
               aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1,
               size = 1.75, label.size = NA, label.padding = unit(0.5, 'lines'),
               fill = NA) +
    geom_segment(aes(x = mz, xend = mz, y = intens, yend = 0, color = match),
                 size = 0.35) +
    geom_hline(aes(yintercept = 0), size = 0.25, color = 'grey20') +
    scale_x_continuous('m/z', expand = c(0.1, 0)) +
    scale_y_continuous('Intensity', limits = c(-1.5, 1.5)) +
    scale_color_manual('', values = pal) +
    boxed_theme(size_sm = 5, size_lg = 6) +
    theme(axis.line.x = element_line(size = 0.3, color = 'grey40'),
          axis.ticks.x = element_line(size = 0.3, color = 'grey40'),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.spacing = unit(1, 'lines'),
          legend.key.width = unit(0.4, 'lines'),
          legend.position = 'none',
          plot.margin = margin(rep(0, 4)),
          plot.title = element_text(size = 6, hjust = 0.5, margin = margin(rep(2, 4))),
          plot.subtitle = element_text(size = 5, hjust = 0.5)) +
    ggh4x::force_panelsizes(rows = unit(2.4, 'cm'),
                            cols = unit(3, 'cm'))
  p1
})
p = ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 4)
p
ggsave("fig/final/figure5b.pdf", p, width = 17 * 4 / 5, height = 8, 
       units = "cm", useDingbats = FALSE)
