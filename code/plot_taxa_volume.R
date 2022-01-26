# code for analysis not included in paper
#requires "processed/psrelmelted.rds" generated from preprocessing.R

library(phyloseq)
library(tidyr)
library(ggplot2)
library(here)

ps.rel.melted <- readRDS(file=here("processed/psrelmelted.rds"))

df <- data.frame(volume = c(), unqs = c())
#number of unique taxa by sample
df <- bind_rows(df, ps.rel.melted %>% group_by(sample) %>%
                  summarise(volume = 0.06*n_distinct(sample), unqs = sum(Abundance > 0), ) %>% select(volume, unqs))
#number of unique taxa by transect
df <- bind_rows(df, ps.rel.melted %>% group_by(site.transect) %>%
                  summarise(volume = 0.06*n_distinct(sample), unqs = sum(Abundance > 0), ) %>% select(volume, unqs))
#number of unique taxa by site
df <- bind_rows(df, ps.rel.melted %>% group_by(site) %>%
                  summarise(volume = 0.06*n_distinct(sample), unqs = sum(Abundance > 0), ) %>% select(volume, unqs))
#number of unique taxa by reef_system
df <- bind_rows(df, ps.rel.melted %>% group_by(reef_system) %>%
                  summarise(volume = 0.06*n_distinct(sample), unqs = sum(Abundance > 0), ) %>% select(volume, unqs))

p <- ggplot(df, aes(x=volume, y=unqs)) +
  geom_point() + labs(y="Unique ASVs", x="Volume seawater sampled (L)") +
  geom_smooth(method="lm", formula=y~x) +
  theme_light()

ggsave(p, filename = here("figures/taxa_volume.pdf"), width = 5, height = 5)
