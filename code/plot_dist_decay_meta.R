# code for plotting the distance decay graph of the secondary analysis Figure 3C
# also analyzes the relationship between temperature/depth and dissimilarity
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)
library(tidyverse)
library(here)

source(here("code/plotting_helpers.R"))

taxonomy <- read.delim(here("SMeta/processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("SMeta/processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("SMeta/processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

# making the imputed abundance counts
f <- zCompositions::cmultRepl(t(asvs), method="CZM", label=0) %>% t()
ps <- phyloseq(otu_table(f, taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
# clr transform
ps.clr <- ps %>% microbiome::transform("clr")
# aitchison distance generation
ps.clr.aitch <- distance(ps.clr, method="euclidean")
# geographic distance matrix
geo <- microbiome::meta(ps) %>% select(lat, lon) %>% rownames_to_column(var="name")
geo.dist <- round(GeoDistanceInMetresMatrix(geo)/1000)
# temperature and depth distance matrix
temp.dist <- stats::dist(metadata %>% select(temperature) %>% filter(!is.na(temperature))) %>% as.matrix()
depth.dist <- stats::dist(metadata %>% select(collectionDepth) %>% filter(!is.na(collectionDepth))) %>% as.matrix()
# pairwise distances for community and geographic/temperature/depth distance
decay.site <- ps.clr.aitch %>% as.matrix() %>% as.data.frame.table(responseName="aitch_distance") %>%
  left_join(as.data.frame.table(geo.dist, responseName="geo_distance")) %>%
  left_join(as.data.frame.table(temp.dist, responseName="temp_distance")) %>%
  left_join(as.data.frame.table(depth.dist, responseName="depth_distance"))
# linear fit line for geographic distance vs aitchison distance
m_lm <- summary(lm(formula = geo_distance ~ aitch_distance, data = decay.site))
eqs <- lm_eqn(m_lm) # getting a printable equation

# Plotting the distance decay
p <- ggplot(decay.site, aes(x=geo_distance, y=-aitch_distance)) + geom_point()
p <- p + theme_light() + scale_x_log10() +
  labs(x="Geographic distance (km)", y="Similarity (-Aitchison distance)") +
  geom_smooth(method='lm', formula=y~x, se=FALSE) +
  theme(axis.title = element_text(size=8), plot.title = element_text(size=8, hjust = .5))
p +
  annotate("label", x=100, y=seq(-120,-155,length=3), label=eqs, parse=TRUE, size=4)

ggsave(filename=here("figures/dist_decay_meta.png"), width = 81, height=81, units = "mm")
# geographic vs aitchison distance mantel correlation
m_mantel <- mantel(ps.clr.aitch, geo.dist)
# Mantel tests require no NA values, so these pipes remove the samples for which temp or depth data are unavailable
temp_mantel <- mantel(ps.clr %>% subset_samples(!is.na(temperature)) %>% distance("euclidean"),
                      temp.dist)
depth_mantel <- mantel(ps.clr %>% subset_samples(!is.na(collectionDepth)) %>% distance("euclidean"),
                       depth.dist)

# saving the mantel test and lm results
saveRDS(m_mantel, file=here("processed/m_mantel.rds"))
saveRDS(temp_mantel, file=here("processed/temp_mantel.rds"))
saveRDS(depth_mantel, file=here("processed/depth_mantel.rds"))
saveRDS(m_lm, file=here("processed/m_lm.rds"))
