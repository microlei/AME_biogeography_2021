# code for plotting figure 4, distance decay measurements at multiple spatial scales
library(ggplot2)
library(cowplot)
library(here)

source(here("code/plot_dist_decay_spatial.R"))
source(here("code/plot_dist_decay_meta.R"))
source(here("code/plot_dist_decay_transect.R"))

a <- ggdraw() + draw_image(here("figures/dist_decay_transect.png"))
b <- ggdraw() + draw_image(here("figures/dist_decay_spatial.png"))
c <- ggdraw() + draw_image(here("figures/dist_decay_meta.png"))

plot_grid(a, b, c, labels=c("A", "B", "C"), nrow=1)
ggsave(filename=here("figures/dist_decay_combined.pdf"),width=169, height=85, units="mm")
