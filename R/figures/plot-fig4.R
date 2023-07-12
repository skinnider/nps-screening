# Plot Figure 4.
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

# plot an example over time
samples_per_month = dates %>% 
  mutate(month = lubridate::floor_date(date, "month")) %>% 
  dplyr::count(month, name = 'total')
compounds_100 = names(which(table(best2$Compound.Name) >= 100))
compounds_20 = table(best2$Compound.Name) %>% 
  sort() %>% 
  extract(. >= 20) %>% 
  names()
time = best2 %>% 
  filter(Compound.Name %in% compounds_20) %>% 
  mutate(month = lubridate::floor_date(date, "month")) %>% 
  group_by(Compound.Name, month) %>% 
  summarise(n = sum(matches > 0)) %>% 
  ungroup() %>% 
  left_join(samples_per_month, by = 'month') %>% 
  mutate(pct = n / total) %>% 
  tidyr::complete(Compound.Name, month, fill = list(n = 0, pct = 0))
pal = c(pals::stepped()[seq(1, 20, 4)],
        pals::stepped3()[seq(1, 20, 4)])
p = time %>% 
  mutate(Compound.Name = factor(Compound.Name, levels = rev(compounds_20))) %>% 
  ggplot(aes(x = month, y = pct, fill = Compound.Name, color = Compound.Name)) +
  facet_wrap(~ Compound.Name, scales = 'free_y', nrow = 2) +
  geom_point(size = 0.8, shape = 21, color = 'black', stroke = 0.1) + 
  geom_smooth(method = 'loess', alpha = 0.15, size = 0.3) +
  scale_x_date('Date', limits = range(dates$date)) +
  scale_y_continuous('% of samples', labels = ~ . * 100,
                     limits = c(-0.001, NA)) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  boxed_theme() +
  theme(legend.position = 'none',
        aspect.ratio = 0.8)
p
ggsave("fig/final/figure4.pdf", p,
       width = 16, height = 6, units = "cm", useDingbats = FALSE)
