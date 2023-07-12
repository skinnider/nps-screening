setwd("~/git/nps-screening")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(magrittr)
source("R/theme.R")

# read performance
perf = readRDS("data/highresnps/performance.rds")

###############################################################################-
## Overall performance ####
###############################################################################-

long = perf %>% 
  pivot_longer(sensitivity:F1, names_to = 'metric') %>% 
  filter(metric != 'F1') %>% 
  filter(is.finite(value))
means = long %>% 
  group_by(metric) %>% 
  summarise(mean = mean(value), median = median(value)) %>% 
  ungroup()
xlabs = with(means, setNames(paste0(str_to_title(metric), '\n(', 
                                    format(median, digits = 2), ')'), metric))
pal = grafify::graf_palettes$contrast %>% unname %>% extract(c(4, 5, 1))
p1 = long %>% 
  mutate(metric = fct_relevel(metric, 'sensitivity', 'specificity', 'accuracy')) %>% 
  ggplot(aes(x = metric, y = value, color = metric, fill = metric)) +
  geom_boxplot(alpha = 0.4, width = 0.6, outlier.shape = NA, size = 0.35) +
  geom_jitter(shape = 21, height = 0, width = 0.15, stroke = 0.15, size = 0.4,
              fill = NA) +
  scale_x_discrete(labels = xlabs) +
  scale_y_continuous('Value', limits = c(0, 1)) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  boxed_theme() +
  theme(aspect.ratio = 1.8,
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p1

###############################################################################-
## By drug class ####
###############################################################################-

# plot performance by drug class
table(long$DrugClass) %>% sort
long %<>%
  mutate(class = ifelse(DrugClass %in% c('Piperidines & pyrrolidines',
                                         'Plants & extracts',
                                         'Indolalkylamines',
                                         'Unknown'),
                        'Other', DrugClass))
means2 = long %>% 
  group_by(class, metric) %>% 
  summarise(mean = mean(value)) %>% 
  ungroup()
pal2 = grafify::graf_palettes$fishy %>% unname
pal2 = pals::kelly()[c(6, 7, 8, 10, 9, 11)]
p2 = long %>% 
  filter(metric != 'accuracy') %>% 
  ggplot(aes(x = class, y = value, color = class, fill = class)) +
  facet_grid(metric ~ ., scales = 'free_x', labeller = as_labeller(str_to_title)) +
  geom_boxplot(alpha = 0.4, width = 0.6, outlier.shape = NA, size = 0.35) +
  geom_jitter(shape = 21, height = 0, width = 0.15, stroke = 0.15, size = 0.4,
              fill = NA) +
  scale_x_reordered() +
  scale_y_continuous('Value', limits = c(0, 1)) +
  scale_fill_manual(values = pal2) +
  scale_color_manual(values = pal2) +
  boxed_theme() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p2

###############################################################################-
## Combine and assemble ####
###############################################################################-

p = p1 | p2
p
ggsave("fig/final/figure6.pdf", p, width = 8, height = 8,
       units = "cm", useDingbats = FALSE)
