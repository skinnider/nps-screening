# Plot Figure 3.
setwd("~/git/nps-screening")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(Spectra)
source("R/theme.R")

# read data
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

# calculate level of evidence for each peak
fragments = databases$NPS %>% 
  filter(MS.Order == 'ms2') %>% 
  distinct(Compound.Name, Product.m.z)
evidence_per_fragment = ms2 %>% 
  dplyr::select(-Product.m.z) %>% 
  left_join(fragments, by = 'Compound.Name') %>% 
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

# tag ppm
ms1 %<>%
  mutate(ppm_err = 1e6 * (mz - m.z) / m.z)

# keep only best evidence per file
best = ms1 %>% 
  group_by(file, Compound.Name) %>% 
  arrange(desc(matches)) %>% 
  dplyr::slice(1) %>% 
  ungroup()

# filter to 2+ fragments
best2 = filter(best, matches >= 2)

# pick spectra for 10 compounds
manual_spectra = c(
  'Eutylone' = 'E3340075726',
  'Fluorofentanyl' = 'R3530202405',
  'Alpha-PiHP' = 'E3200380919',
  '8-Aminoclonazolam' = 'R3400205756',
  'Bromazolam' = 'R3550029494',
  'N-Ethylpentedrone' = 'R3410200601',
  'Napthyl-U-47700' = 'R3480306158',
  'Metizolam' = 'R3370230148',
  'Furanyl UF-17' = 'R3380109563',
  'N-Ethyl-U-47700' = 'E3200443744'
) %>%
  data.frame(Compound.Name = names(.), file = .)
sample = inner_join(ms1, manual_spectra) %>% 
  filter(spectrum %in% best2$spectrum)
distinct(sample, file, Compound.Name, spectrum) 

# plot MS/MS for the sample
msms = ms2 %>% 
  dplyr::select(-(Compound.Name:Integration.Strategy)) %>% 
  inner_join(sample, by = c('file', 'spectrum', 'date'))
plots = map(unique(sample$Compound.Name), ~ {
  compound = .x
  msms0 = filter(msms, Compound.Name == compound) %>% 
    # rename fragment ion m/z
    dplyr::rename(mz = mz.x) %>% 
    # sqrt-transform
    mutate(i = sqrt(i)) %>%
    # normalize fragment intensities
    group_by(file, spectrum, date) %>% 
    mutate(intens = i / max(i)) %>%
    # sqrt-transform
    # mutate(intens = sqrt(intens)) %>% 
    ungroup()
  
  # get reference
  ref = databases$NPS %>% 
    filter(Workflow == 'Fragment', Compound.Name == compound) %>% 
    dplyr::rename(mz = Product.m.z) %>% 
    mutate(intens = -0.5)
  
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
    ggtitle(compound) +
    geom_blank() + 
    geom_rect(data = bg,
              aes(xmin = -Inf, xmax = Inf, ymax = 0, ymin = -Inf),
              alpha = 0.1, fill = 'grey70', color = NA) +
    geom_label(data = bg %>% mutate(label = 'Database'),
               aes(x = Inf, y = -Inf, label = label), hjust = 1, vjust = 0,
               size = 1.75, label.size = NA, label.padding = unit(0.5, 'lines'),
               fill = NA) +
    geom_label(data = bg %>% mutate(label = 'UDS'),
               aes(x = -Inf, y = Inf, label = label), hjust = 0, vjust = 1,
               size = 1.75, label.size = NA, label.padding = unit(0.5, 'lines'),
               fill = NA) +
    geom_segment(aes(x = mz, xend = mz, y = intens, yend = 0, color = match),
                 size = 0.35) +
    geom_text(data = filter(mirror_df, intens < 0),
              aes(x = mz, y = intens, color = match, # hjust = -sign(intens), 
                  label = formatC(mz, format = 'f', digits = 2) %>% trimws),
              size = 1.5, angle = 45, hjust = 1, nudge_y = -0.1) +
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
          plot.margin = margin(rep(0, 4))) +
    ggh4x::force_panelsizes(rows = unit(2.4, 'cm'),
                            cols = unit(3, 'cm'))
  p1
})
p = ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 5)
p
ggsave("fig/final/figure3.pdf", p, width = 17, height = 8, 
       units = "cm", useDingbats = FALSE)
