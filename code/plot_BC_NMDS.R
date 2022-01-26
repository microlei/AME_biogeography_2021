# Code for plotting the Bray-Curtis NMDS and beta-dispersion for Figure 1
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggforce)
library(cowplot)
library(vegan)
library(here)
source(here("code/plotting_helpers.R"))

# loading files
taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
sample_data(ps)$site <- factor(sample_data(ps)$site, levels = breaks.site) # make factor labels pretty
ps.rel <- ps %>% transform_sample_counts(function(x) 100*x/sum(x))

#BC NMDS
dist.rel <- phyloseq::distance(ps.rel, method="bray")
nmds <- metaMDS(dist.rel, trymax=100)
nmds_plot <- plot_ordination(ps.rel, nmds, color="site", shape="site")
nmds_plot <- nmds_plot + geom_point(aes(color=site, shape=site), size=4) + scale_color_manual(values=cols.site, labels=labels.site, breaks=breaks.site, name="Reef") +
  scale_shape_manual(values=c(15, 16, 17, 3, 4, 9), labels=labels.site, breaks=breaks.site, name="Reef") +
  theme_light() + theme(legend.text=element_text(size=8)) +
  labs(xlab="NMDS1", ylab="NMDS2")
nmds_legend <- get_legend(nmds_plot + guides(color = guide_legend(nrow=2))+ theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))
nmds_plot <- nmds_plot+theme(legend.position = "none")
nmds_plot
# beta dispersion
site <- sample_data(ps.rel)$site
disp.rel.site <- betadisper(d = dist.rel, group=site)
df <- data.frame(disper=disp.rel.site$distances,group=disp.rel.site$group)
disper_plot<-ggplot(df, aes(x=group, y=disper)) + geom_boxplot(aes(color=group)) + scale_color_manual(values=cols.site, labels=labels.site, breaks=breaks.site, name=NULL) +
  labs(y="Distance from centroid") +coord_flip() + scale_x_discrete(limits=rev) +
  geom_text(data=tibble(x=c(1,2,3,4,5,6), y=rep(.6,6), sig=c("a","b","ab","ab","b","ab")), aes(x=x,y=y,label=rev(sig)), size=5) + # significance groups from disper_adonis_tests.R
  theme_light() + theme(axis.title.y.left = element_blank(),
                        axis.text.y.left = element_blank(),
                        legend.position = "none")
disper_plot
# original layout
combined <- plot_grid(nmds_plot, disper_plot, labels = "AUTO") %>% plot_grid(nmds_legend, ncol=1, rel_heights = c(1,.1))
ggsave(combined, filename = here("figures/BC_NMDS.pdf"), width=160, height=120, units="mm")

# different layout to test <- ended up going with this one for publication
nmds_legend <- get_legend(disper_plot + guides(color=guide_legend(ncol=1, byrow= T)) + theme(legend.position = "left", legend.box.margin = margin(0,0,0,0,"mm"), legend.spacing = unit(0,"cm"), legend.title = element_blank(), legend.key.height = unit(1.5,"cm"), legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(2, "mm")))
plot_grid(nmds_plot, nmds_legend, disper_plot, labels=c("A", "B", ""), nrow=1, rel_widths = c(1,.6,1))
ggsave(filename = here("figures/BC_NMDS_test.pdf"), width=160, height=120, units="mm")
