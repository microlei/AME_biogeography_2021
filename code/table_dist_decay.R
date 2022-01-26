# code for making the table of mantel and slope values
library(vegan)
library(tidyverse)
library(here)

#depends on plot_dist_decay_spatial.R and plot_dist_decay_meta.R

sp_mantel <- readRDS(file=here("processed/sp_mantel.rds"))
m_mantel <- readRDS(file=here("processed/m_mantel.rds"))
fla_mantel <- readRDS(file=here("processed/fla_mantel.rds"))
stt_mantel <- readRDS(file=here("processed/stt_mantel.rds"))

sp_lm <- readRDS(file=here("processed/sp_lm.rds"))
m_lm <- readRDS(file=here("processed/m_lm.rds"))
fla_lm <- readRDS(file=here("processed/fla_lm.rds"))
stt_lm <- readRDS(file=here("processed/stt_lm.rds"))
trans_lm <- readRDS(file=here("processed/trans_lm.rds"))

t<-tibble(
  Scope = c("Within transects", "(Reef system) USVI", "(Reef system) Florida Keys", "All samples in Fl/VI-based study", "All samples in secondary analysis"),
  Scale = c("0-9 m", "0-3 km", "0-279 km", "0-1,978 km", "0-16,874 km"),
  Mantel.r = format(c(NA, stt_mantel$statistic, fla_mantel$statistic, sp_mantel$statistic, m_mantel$statistic),digits=4, scientific=F),
  Mantel.pvalue = c(NA, stt_mantel$signif, fla_mantel$signif, sp_mantel$signif, m_mantel$signif),
  Slope = format(c(-trans_lm$coefficients[2], -stt_lm$coefficients[2], -fla_lm$coefficients[2], -sp_lm$coefficients[2], -m_lm$coefficients[2]), digits=1, scientific = F),
)

colnames(t) <- c("Scope", "Spatial scale", "Correlation (Mantel r)", "Mantel p value", "Slope of linear fit")
readr::write_delim(t, file=here("figures/table_dist_decay.csv"), delim=",")
