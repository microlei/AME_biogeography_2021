# code for plotting figure 1B
library(ggplot2)
library(here)
source(here("code/plotting_helpers.R"))

world <- map_data("world")
worldmeta <- read.delim(here("SMeta/output/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill="grey", alpha=0.3) +
  geom_point(data=worldmeta, aes(x=lon, y=lat, color=study, shape=study), size=3) + scale_color_manual(values=cols.study, labels=labels.study) +
  scale_shape_manual(values = c(0,1,2,5,6), labels=labels.study) + coord_map(ylim = c(-25,35)) +
  theme_light() + theme(legend.position = "bottom", legend.title = element_blank(), legend.margin=margin(0,0,0,0), plot.margin=margin(0,0,0,0)) +
  labs(x="", y="") +
  guides(color=guide_legend(nrow=2, byrow=TRUE))

ggsave(filename = here("figures/map_meta.png"), width = 169, height = 60, units="mm")
