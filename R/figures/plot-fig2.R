# Plot Figure 3.
setwd("~/git/nps-screening")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
library(Spectra)
library(MSnbase)
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

###############################################################################-
# Bar chart: number of fragments ####
###############################################################################-

pal1 = c('grey88', pals::stepped3()[c(4, 3, 2, 1) + 0]) %>% setNames(1:5)
n = dplyr::count(best2, Compound.Name)
p1 = best2 %>% 
  ggplot(aes(x = fct_infreq(Compound.Name), fill = factor(matches))) +
  geom_bar(color = 'grey10', size = 0.15, width = 0.7) + 
  geom_text(data = n, aes(y = n, label = n, fill = NA), size = 1.5, vjust = 0,
            nudge_y = 20) +
  scale_fill_manual('# of matching product ions', values = pal1) +
  scale_y_continuous('# of samples', expand = expansion(c(0, 0.1)),
                     labels = ~ replace(., . == 500, '>500')) +
  coord_cartesian(ylim = c(0, 500)) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, 'lines'),
        aspect.ratio = 0.3)
p1

###############################################################################-
# Bar chart: proportion of fragments ####
###############################################################################-

pal2 = pals::cubehelix(100) %>% setNames(seq(1, 100))
# pal = pals::kovesi.rainbow(100) %>% setNames(seq(1, 100))
p2 = best2 %>% 
  mutate(pct = round(100 * matches / total_fragments)) %>% 
  ggplot(aes(x = fct_infreq(Compound.Name), fill = factor(pct))) +
  geom_bar(color = 'grey10', size = 0.15, width = 0.7) + 
  scale_fill_manual('% of matching\nproduct ions', values = pal2,
                    limits = seq_len(100),
                    breaks = seq(10, 100, 10),
                    labels = seq(10, 100, 10)) +
  scale_y_continuous('# of samples', expand = expansion(c(0, 0.1)),
                     labels = ~ replace(., . == 500, '>500')) +
  coord_cartesian(ylim = c(0, 500)) +
  # guides(fill)
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.35, 'lines'),
        legend.position = 'right',
        aspect.ratio = 0.3)
p2

###############################################################################-
# Boxplot: retention times ####
###############################################################################-

# retention times
ranges = databases$NPS %>% 
  filter(Compound.Name %in% best2$Compound.Name,
         MS.Order == 'ms1') %>% 
  mutate(Compound.Name = factor(Compound.Name, 
                                levels = levels(fct_infreq(best2$Compound.Name))))
p3 = best2 %>% 
  ggplot(aes(x = fct_infreq(Compound.Name), y = rt)) +
  geom_blank() +
  geom_rect(data = ranges, 
            aes(xmin = as.integer(fct_infreq(Compound.Name)) - 0.4, 
                xmax = as.integer(fct_infreq(Compound.Name)) + 0.4, 
                ymin = 60 * Retention.Time - Retention.Time.Window,
                ymax = 60 * Retention.Time + Retention.Time.Window,
                y = 60 * Retention.Time),
            alpha = 0.8, fill = pals::stepped()[4], color = NA, # 'black', 
            size = 0.2) +
  geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA,
               fill = 'grey82') +
  geom_jitter(height = 0, width = 0.3, shape = 21, size = 0.4, stroke = 0.25) +
  scale_y_continuous('Retention time (s)', limits = range(ms1$rt)) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, 'lines'),
        aspect.ratio = 0.3)
p3

###############################################################################-
# Boxplot: mass accuracy ####
###############################################################################-

p4 = best2 %>% 
  ggplot(aes(x = fct_infreq(Compound.Name), y = ppm_err)) +
  geom_hline(aes(yintercept = 0), color = 'grey88', size = 0.45) +
  geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA,
               fill = 'grey82') +
  geom_jitter(height = 0, width = 0.3, shape = 21, size = 0.4, stroke = 0.25) +
  scale_y_continuous('Mass accuracy (ppm)', breaks = seq(-8, 10, 4)) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, 'lines'),
        aspect.ratio = 0.4)
p4

###############################################################################-
# Boxplot: MS/MS similarity ####
###############################################################################-

## sample 100 spectra per compound
set.seed(0)
sample = best2 %>% 
  group_by(Compound.Name) %>% 
  filter(n() > 1) %>% 
  sample_n(min(100, n())) %>% 
  ungroup() %>% 
  distinct()

# calculate dot-product
dot_product_file = "data/figures/dot-products.rds"
if (file.exists(dot_product_file)) {
  dot_products = readRDS(dot_product_file)
} else {
  dot_products = map_dfr(unique(sample$Compound.Name), ~ {
    compound = .x
    message(compound)
    ## randomly sample a maximum of 100
    set.seed(0)
    sample = best2 %>% 
      filter(Compound.Name == compound) %>% 
      distinct(file, spectrum, date, rt) %>% 
      sample_n(min(100, n())) %>% 
      distinct()
    ## create spectra
    msms = ms2 %>% 
      dplyr::select(-(Workflow:Integration.Strategy)) %>% 
      inner_join(sample, by = c('file', 'spectrum', 'date'))
    msms0 = filter(msms, Compound.Name == compound) %>% 
      # rename fragment ion m/z
      dplyr::rename(mz = mz) %>% 
      # normalize fragment intensities
      group_by(file, spectrum, date) %>% 
      mutate(intens = i / max(i)) %>% 
      ungroup()
    
    ## create spectra
    spd = distinct(msms0, spectrum, rt) %>% 
      mutate(msLevel = 2L) %>% 
      dplyr::rename(rtime = rt) %>% 
      arrange(spectrum)
    spd$mz = msms0 %>% split(.$spectrum) %>% map('mz')
    spd$intensity = msms0 %>% split(.$spectrum) %>% map('intens')
    stopifnot(all(spd$spectrum == names(spd$mz))) ## TRUE
    stopifnot(all(spd$spectrum == names(spd$intensity))) ## TRUE
    sps = Spectra(spd)
    sps$name = spd$spectrum
    dp = compareSpectra(sps, ppm = 20)
    reshape2::melt(dp, varnames = c('spectrum1', 'spectrum2'),
                   value.name = 'dot_product') %>% 
      filter(spectrum1 < spectrum2) %>% 
      mutate(compound = compound)
  })
  saveRDS(dot_products, dot_product_file)
}
# plot
levels = fct_infreq(best2$Compound.Name) %>% 
  levels()
p5 = dot_products %>% 
  ggplot(aes(x = factor(compound, levels = levels), y = dot_product)) +
  geom_boxplot(size = 0.35, alpha = 0.4, width = 0.6, outlier.shape = NA,
               fill = 'grey82') +
  geom_jitter(height = 0, width = 0.3, shape = 21, size = 0.3, stroke = 0.1) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous('Dot product', limits = c(0, 1)) +
  boxed_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.key.size = unit(0.4, 'lines'),
        aspect.ratio = 0.3)
p5

###############################################################################-
# Heatmap: fragment overlap ####
###############################################################################-

# create fragment matrices for suspect database
fragments = databases$NPS %>% 
  filter(MS.Order == 'ms2') %>% 
  distinct(Compound.Name, Product.m.z)
nps = unique(fragments$Compound.Name)
matrices = map(nps, ~ {
  mat = filter(fragments, Compound.Name == .x) %>% 
    arrange(Product.m.z) %>% 
    distinct() %>% 
    pull(Product.m.z) %>% 
    cbind(1)
}) %>% setNames(nps)

# read HighResNPS
highresnps = read_excel("data/databases/highresnps january 2023 consensus.xlsx")
nps_highresnps = unique(highresnps$Compound) %>% setdiff('')
# create fragment matrices
fragments_highresnps = highresnps %>% 
  dplyr::select(Compound, ends_with('mass')) %>% 
  pivot_longer(ends_with('mass')) %>% 
  drop_na(value)
matrices_highresnps = map(nps_highresnps, ~ {
  mat = filter(fragments_highresnps, Compound == .x) %>% 
    arrange(value) %>% 
    distinct() %>% 
    pull(value) %>% 
    cbind(1)
}) %>% setNames(nps_highresnps)

pairs_file = "data/spectral_similarity/suspect-vs-highresnps.rds"
if (!file.exists(pairs_file2)) {
  pairs_highresnps = tidyr::crossing(nps_x = nps, nps_y = nps_highresnps) %>% 
    filter(nps_x != nps_y) %>% 
    pmap_dfr(function(...) {
      current = tibble(...)
      message(paste(current, collapse = ' / '))
      mat_x = matrices[[current$nps_x]]
      mat_y = matrices_highresnps[[current$nps_y]]
      n_overlap = joinPeaks(mat_x, mat_y, ppm = 20, type = "inner")$x %>% 
        nrow()
      mutate(current, 
             overlap = n_overlap,
             n_x = nrow(mat_x),
             n_y = nrow(mat_y))
    })
  saveRDS(pairs_highresnps, pairs_file)
} else {
  pairs_highresnps = readRDS(pairs_file)
}

# nearest-neighbors
nn_highresnps = pairs_highresnps %>% 
  # remove overlaps
  filter(nps_x %in% best2$Compound.Name) %>% 
  filter(!(nps_x == '2-Fluorodeschloroketamine' & nps_y == '2-Fluoro Deschloroketamine'),
         !(nps_x == '4-Trifluoromethyl-U-47700' & nps_y == '4-TFM U-47700'),
         !(nps_x == 'Cl-PCP' & nps_y == '3-Cl-PCP'),
         !(nps_x == 'F-PCP' & nps_y == '3F-PCP'),
         !(nps_x == 'Fluorofentanyl' & nps_y == 'meta-Fluoro Fentanyl'),
         !(nps_x == 'MeO-PCE' & nps_y == '3-MeO-PCE'),
         !(nps_x == 'N-methyl-U-47931E' & nps_y == 'N-Methyl U-47931E'),
         !(nps_x == 'n-Piperidinyl etonitazene' & nps_y == 'N-Piperidinyl Etonitazene'),
         !(nps_x == 'Napthyl-U-47700' & nps_y == '1-Naphthyl U-47700'),
         !(nps_x == 'Napthyl-U-47700' & nps_y == '2-Naphthyl U-47700')
  ) %>% 
  mutate(pct_overlap = overlap / n_x) %>% 
  group_by(nps_x) %>% 
  arrange(desc(pct_overlap), desc(overlap)) %>% 
  dplyr::slice(1) %>% 
  ungroup() 

# prepare for plotting
df6 = nn_highresnps %>% 
  # reorder to match detection frequency
  filter(nps_x %in% best2$Compound.Name) %>% 
  mutate(nps_x = factor(nps_x, levels = levels(fct_infreq(best2$Compound.Name))))

# plot
pal_num = c('white', 'white', pals::stepped3()[c(4, 3, 2, 1) + 0]) %>% setNames(c(0:4, 6))
p6 = df6 %>% 
  ggplot(aes(x = nps_x, y = 'Fragment\noverlap', fill = factor(overlap))) +
  geom_tile(color = 'white') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual('# of matching product ions', values = pal_num,
                    breaks = c('1', '2', '3', '4', '6'),
                    labels = c('0-1', '2', '3', '4', '5+')) +
  guides(fill = guide_legend(override.aes = list(color = 'black'))) +
  coord_fixed() + 
  boxed_theme() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.35, 'lines'))
p6

###############################################################################-
# Combine and assemble ####
###############################################################################-

# revised version
p = (p1 + theme(axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(),
                aspect.ratio = 0.25)) +
  (p2 + theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              aspect.ratio = 0.25)) +
  (p6 + theme(axis.text.x = element_blank(),
               legend.position = 'none')) +
  (p3 + theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              aspect.ratio = 0.25)) +
  (p4 + theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              aspect.ratio = 0.25)) +
  (p5 + theme(aspect.ratio = 0.25)) + 
  plot_layout(ncol = 1) &
  theme(plot.margin = margin(rep(0.6, 4)),
        plot.background = element_blank())
p
ggsave("fig/final/figure2-revised.pdf", p, width = 12, height = 13, units = "cm", 
       useDingbats = FALSE)   
