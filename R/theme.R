library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(patchwork)
library(paletteer)
library(scico)
library(drlib)
library(BuenColors)
library(ggrastr)

# Define theme
clean_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
        axis.text.y = element_text(size = size_sm),
        axis.ticks.length.x = unit(0.15, 'lines'),
        axis.ticks.length.y = unit(0.15, 'lines'),
        axis.title.x = element_text(size = size_lg),
        axis.title.y = element_text(size = size_lg),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size = size_sm),
        # strip.background = element_rect(fill = "grey90", color = "grey90",
        #                                 size = 0),
        strip.background = element_blank(),
        axis.line.y = element_line(colour = "grey50"),
        axis.line.x = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.position = "top",
        legend.text = element_text(size = size_sm),
        legend.title = element_text(size = size_sm),
        legend.key.size = unit(0.6, "lines"),
        legend.margin = margin(rep(0, 4)),
        # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
        # legend.box.spacing = ggplot2::margin(rep(0, 4)),
        legend.background = element_blank(),
        plot.title = element_text(size = size_lg, hjust = 0.5),
        axis.ticks.length = unit(2, 'pt'),)
}

grid_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.ticks.length.x = unit(0.15, 'lines'),
          axis.ticks.length.y = unit(0.15, 'lines'),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          strip.text = element_text(size = size_sm),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_line(colour = "grey50"),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
          # legend.box.spacing = ggplot2::margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5),
          axis.ticks.length = unit(2, 'pt'),)
}

boxed_theme = function(size_lg = 6, size_sm = 5) {
  theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = size_sm),
          axis.text.y = element_text(size = size_sm),
          axis.ticks.length.x = unit(0.15, 'lines'),
          axis.ticks.length.y = unit(0.15, 'lines'),
          axis.title.x = element_text(size = size_lg),
          axis.title.y = element_text(size = size_lg),
          panel.grid = element_blank(),
          strip.text = element_text(size = size_sm),
          strip.background = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks = element_line(colour = "grey50"),
          legend.position = "top",
          legend.text = element_text(size = size_sm),
          legend.title = element_text(size = size_sm),
          legend.key.size = unit(0.6, "lines"),
          legend.margin = margin(rep(0, 4)),
          # legend.box.margin = ggplot2::margin(rep(0, 4), unit = 'lines'),
          # legend.box.spacing = ggplot2::margin(rep(0, 4)),
          legend.background = element_blank(),
          plot.title = element_text(size = size_lg, hjust = 0.5))
}

umap_theme = function(size_lg = 6, size_sm = 5) {
  boxed_theme(size_lg = size_lg, size_sm = size_sm) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_text(hjust = 0, vjust = 0),
          axis.title.x = element_text(hjust = 0))
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}


lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  for (value in col){
    if (value > 255) {
      col[col == value] = 255
    }
  }
  col <- rgb(t(col), maxColorValue=255)
  col
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # remove zero
  l <- gsub("0e\\+00", "0", l)
  # remove one
  l <- gsub("^1e\\+00", "1", l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + from exponent
  l <- gsub("e\\+" ,"e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # remove 1 x 10^ (replace with 10^)
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

plot_pal = function(pal) {
  grid::grid.raster(pal, interpolate=F)
}

cubehelix = function(n_colors) {
  colours = c("#000000", "#1A1935", "#15474E", "#2B6F39", "#767B33", "#C17A6F",
              "#D490C6", "#C3C0F2")
  idxs = 0.3
  if (n_colors > 1)
    idxs = seq(0, 1, 1 / (n_colors - 1))
  colour_ramp(colours)(idxs)
}

kinney6 = c(
  "#c6c3bf",
  "#119e87",
  "#53bad3",
  "#559ed2",
  "#3b5687",
  "#e34d3b")

colours.cafe447 = c('#ffb838', '#fee5a5', '#f7f6fee', '#486d87')
colours.cafe433 = c('#077893', '#e3deca', '#fcfaf1', '#ff9465')
colours.cafe425 = c('#2B5B6C', '#C7CFAC', '#FCFAF1', '#E34F33', '#FFC87E')
colours.cafe322 = c("#7bbaea", "#d46363", "#fbdaa7", "#fcfaf2", "#30598c")

winsorize = function(vec, limits) {
  lower_limit = limits[1]
  upper_limit = limits[2]
  if (!is.na(upper_limit))
    vec[vec > upper_limit] = upper_limit
  if (!is.na(lower_limit))
    vec[vec < lower_limit] = lower_limit
  return(vec)
}

x100 = function(x) x * 100
