# code for plotting figure 3B, distance decay measurements of VI and Fl reefs together and individually
# also performs and saves mantel tests for geographic distance and aitchison distance
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)
library(here)

source(here("code/plotting_helpers.R"))

taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
# imputing zeros
f <- zCompositions::cmultRepl(t(asvs), method="CZM", label=0) %>% t()
ps <- phyloseq(otu_table(as.matrix(f), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
# clr transform
ps.clr <- ps %>% microbiome::transform("clr")

# Distance-decay of all samples in spatial study
ps.aitch <- distance(ps.clr, method="euclidean") # aitchison distance transform
# geographic distance matrix
geo <- microbiome::meta(ps.clr) %>% select(lat, lon) %>% rownames_to_column(var="name")
geo.matrix <- round(GeoDistanceInMetresMatrix(geo) / 1000)
# pairwise distances (geographic and aitchison)
decay.site <- ps.aitch %>% as.matrix() %>% as.data.frame.table(responseName="aitch_distance") %>%
  left_join(as.data.frame.table(geo.matrix, responseName="geo_distance"))
decay.site <- filter(decay.site, aitch_distance!=0)
# lm fit line
sp_lm <- summary(lm(formula = geo_distance ~ aitch_distance, data = decay.site))
eqs <- lm_eqn(sp_lm) # getting variables for lm fit line
# plotting
p <- ggplot(decay.site, aes(x=geo_distance, y=-aitch_distance)) + geom_point()
p <- p + theme_light() + scale_x_log10() +
  labs(x="Geographic distance (km)", y="Similarity (-Aitchison distance)") +
  geom_smooth(method='lm', formula=y~x, se=FALSE) +
  theme(axis.title = element_text(size=8), plot.title = element_text(size=8, hjust = .5))
p +
  annotate("label", x=50, y=seq(-130,-160,length=3), label=eqs, parse=TRUE, size=4)

# saving results for all samples
ggsave(filename=here("figures/dist_decay_spatial.png"), width = 81, height=81, units="mm")
sp_mantel <- mantel(ps.aitch, geo.matrix) # mantel test
saveRDS(sp_mantel, file=here("processed/sp_mantel.rds"))
saveRDS(sp_lm, file=here("processed/sp_lm.rds"))

# Distance-decay by reef system (FLA)
ps.clr.fla <- ps.clr %>% subset_samples(reef_system=="FLA")
ps.aitch.fla <- distance(ps.clr.fla, method="euclidean")
geo.fla <- microbiome::meta(ps.clr.fla) %>% select(lat, lon) %>% rownames_to_column(var="name")
geo.fla <- round(GeoDistanceInMetresMatrix(geo.fla)/1000)
decay.fla <- ps.aitch.fla %>% as.matrix() %>% as.data.frame.table(responseName = "aitch_distance") %>%
  left_join(as.data.frame.table(geo.fla, responseName="geo_distance"))
decay.fla <- filter(decay.fla, aitch_distance!=0)
fla_lm <- summary(lm(formula = geo_distance~aitch_distance, data=decay.fla))
eqs <- lm_eqn(fla_lm)

p <- ggplot(decay.fla, aes(x=geo_distance, y=-aitch_distance)) + geom_point()
p <- p + theme_light() + scale_x_log10() +
  labs(x="Geographic distance (km)", y="Similarity (-Aitchison distance)") +
  geom_smooth(method='lm', formula=y~x, se=FALSE) +
  theme(axis.title = element_text(size=8), plot.title = element_text(size=8, hjust = .5))
p +
  annotate("label", x=50, y=seq(-120,-135,length=3), label=eqs, parse=TRUE, size=2)

#saving results for FLA reef system
ggsave(filename=here("figures/dist_decay_fla.png"), width = 81, height=81, units="mm")
fla_mantel <- mantel(ps.aitch.fla, geo.fla)
saveRDS(fla_mantel, file=here("processed/fla_mantel.rds"))
saveRDS(fla_lm, file=here("processed/fla_lm.rds"))

# Distance-decay by reef system (STT)
ps.clr.stt <- ps.clr %>% subset_samples(reef_system=="STT")
ps.aitch.stt <- distance(ps.clr.stt, method="euclidean")
geo.stt <- microbiome::meta(ps.clr.stt) %>% select(lat, lon) %>% rownames_to_column(var="name")
geo.stt <- round(GeoDistanceInMetresMatrix(geo.stt)/1000)
decay.stt <- ps.aitch.stt %>% as.matrix() %>% as.data.frame.table(responseName = "aitch_distance") %>%
  left_join(as.data.frame.table(geo.stt, responseName="geo_distance"))
decay.stt <- filter(decay.stt, aitch_distance!=0)
stt_lm <- summary(lm(formula = geo_distance~aitch_distance, data=decay.stt))
eqs <- lm_eqn(stt_lm)

p <- ggplot(decay.stt, aes(x=geo_distance, y=-aitch_distance)) + geom_point()
p <- p + theme_light() + scale_x_log10() +
  labs(x="Geographic distance (km)", y="Similarity (-Aitchison distance)") +
  geom_smooth(method='lm', formula=y~x, se=FALSE) +
  theme(axis.title = element_text(size=8), plot.title = element_text(size=8, hjust = .5))
p +
  annotate("label", x=2, y=seq(-120,-135,length=3), label=eqs, parse=TRUE, size=2)

#saving results for STT reef system
ggsave(filename=here("figures/dist_decay_stt.png"), width = 81, height=81, units="mm")
stt_mantel <- mantel(ps.aitch.stt, geo.stt)
saveRDS(stt_mantel, file=here("processed/stt_mantel.rds"))
saveRDS(stt_lm, file=here("processed/stt_lm.rds"))
