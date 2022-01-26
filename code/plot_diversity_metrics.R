# code for figure S2 alpha diversity metrics
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
library(here)

source(here("code/plotting_helpers.R"))

taxonomy <- read.delim("processed/taxonomy.txt", "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim("processed/ASVs.txt", "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim("processed/metadata.csv", ",", header = TRUE, row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
sample_data(ps)$site.transect <- factor(sample_data(ps)$site.transect, levels = breaks.transect)
sample_data(ps)$site <- factor(sample_data(ps)$site, levels=breaks.site)

# plottging the diversity metrics
theme_set(theme_light() + theme(axis.ticks.y=element_blank(),
                                axis.text.y=element_blank(),
                                axis.text.x.bottom=element_text(angle=0, size=8, hjust=.5),
                                axis.title.x.bottom = element_blank(),
                                axis.title.y.left = element_text(size=9),
                                legend.title = element_blank(),
                                legend.text = element_text(size=8),
                                strip.text=element_blank(),
                                legend.position = "none",
                                plot.title=element_text(size=9, hjust=0.5, face = "bold")))

alpha_p <- plot_richness(ps, x="site.transect", measures="Observed", color="site") + geom_boxplot() + scale_x_discrete(limits=rev) + labs(title="Transect (within reef)", x="Observed ASVs") + coord_flip() +
  scale_color_manual(values=cols.site, labels=labels.site, breaks=breaks.site, name="Reef") +
  geom_text(data=tibble(x=15.5, y=1600, sig="*"), aes(x=x, y=y, label=rev(sig)), size=5, inherit.aes = F)
alpha_q <- plot_richness(ps, x="site", measures="Observed", color="site") + geom_boxplot() + scale_x_discrete(limits=rev) + coord_flip() +
  scale_color_manual(values=cols.site) + labs(title="Individual Reef", x="")
alpha_r <- plot_richness(ps, x="reef_system", measures="Observed") + geom_boxplot(aes(fill=reef_system)) + scale_x_discrete(limits=rev) + coord_flip() +
  scale_fill_manual(values=c(cols.site[7], cols.site[8]), labels=c("Florida", "Virgin Islands")) + labs(title="Reef system", x="") +
  geom_text(data=tibble(x=2, y=1600, sig="*"), aes(x=x,y=y,label=sig), size=5, inherit.aes = F)
shannon_p <- plot_richness(ps, x="site.transect", measures="shannon", color="site") + geom_boxplot() + scale_x_discrete(limits=rev) + coord_flip() +
  scale_color_manual(values=cols.site) +labs(x="Shannon Index")
shannon_q <- plot_richness(ps, x="site", measures="shannon", color="site") + geom_boxplot() + scale_x_discrete(limits=rev) + coord_flip() +
  scale_color_manual(values=cols.site)+labs(x="") +
  geom_text(data=tibble(x=c(1,2,3,4,5,6), y=rep(6,6), sig=c("b","a","ab","a","a","b")), aes(x=x,y=y,label=rev(sig)),size=4, inherit.aes = F)
shannon_r <- plot_richness(ps, x="reef_system", measures="shannon") + geom_boxplot(aes(fill=reef_system)) + scale_x_discrete(limits=rev) + coord_flip() +
  scale_fill_manual(values=c(cols.site[7], cols.site[8]))+labs(x="") +
  geom_text(data=tibble(x=2, y=6, sig="*"), aes(x=x,y=y,label=sig), size=5, inherit.aes = F)
simpson_p <- plot_richness(ps, x="site.transect", measures="simpson", color="site") + geom_boxplot() + scale_x_discrete(limits=rev) + coord_flip() +
  scale_color_manual(values=cols.site)+labs(x="Simpson's Index")
simpson_q <- plot_richness(ps, x="site", measures="simpson", color="site") + geom_boxplot() + scale_x_discrete(limits=rev) + coord_flip() +
  scale_color_manual(values=cols.site) +labs(x="") +
  geom_text(data=tibble(x=c(1,2,3,4,5,6), y=rep(1,6), sig=c("c","a","a","a","a","b")), aes(x=x, y=y, label=rev(sig)), size=4, inherit.aes = F)
simpson_r <- plot_richness(ps, x="reef_system", measures="simpson") + geom_boxplot(aes(fill=reef_system)) + scale_x_discrete(limits=rev) + coord_flip() +
  scale_fill_manual(values=c(cols.site[7], cols.site[8]))+labs(x="") +
  geom_text(data=tibble(x=2, y=1, sig="*"), aes(x=x,y=y,label=sig), size=5, inherit.aes = F)

p_legend <- get_legend(alpha_p + guides(color = guide_legend(nrow = 1))+ theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))
r_legend <- get_legend(alpha_r + guides(color = guide_legend(nrow = 1))+ theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))
patchwork <- (alpha_p + alpha_q + alpha_r)/(shannon_p + shannon_q + shannon_r) /(simpson_p + simpson_q + simpson_r) /p_legend /r_legend
patchwork <- patchwork + plot_layout(heights = c(1,1,1,.1,.1))
ggsave(filename=here("figures/diversity.pdf"), plot=patchwork, width=160, height=150, unit="mm")
