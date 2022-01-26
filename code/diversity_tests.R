## Code for testing the significance of differences in alpha diversity of transect-based study
library(phyloseq)
library(vegan)
library(tidyverse)
library(here)
source(here("code/plotting_helpers.R"))

# load files
taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
sample_data(ps)$site.transect <- factor(sample_data(ps)$site.transect, levels = breaks.transect)
sample_data(ps)$site <- factor(sample_data(ps)$site, levels=breaks.site)

# significance values for alpha diversity metrics at multiple spatial scales
rich_est <- estimate_richness(ps, measures=c("Observed", "Shannon", "Simpson"))
rich_est <- rich_est %>% rownames_to_column(var="sample") %>% left_join(metadata, by="sample") %>% select(sample, Observed, Shannon, Simpson, site, reef_system, transect, site.transect)

# pairwise t-tests for alpha diversity metrics
# these results are then used in figure 2 to make the significance groups
transect.simpson <- rich_est %>% select(site, site.transect, Simpson) %>% group_by(site) %>% summarise(t_test=list(pairwise.t.test(Simpson, site.transect)))
transect.shannon <- rich_est %>% select(site, site.transect, Shannon) %>% group_by(site) %>% summarise(t_test=list(pairwise.t.test(Shannon, site.transect)))
transect.observed <- rich_est %>% select(site, site.transect, Observed) %>% group_by(site) %>% summarise(t_test=list(pairwise.t.test(Observed, site.transect)))

site.simpson <- pairwise.t.test(rich_est$Simpson, rich_est$site)
site.shannon <- pairwise.t.test(rich_est$Shannon, rich_est$site)
site.observed <- pairwise.t.test(rich_est$Observed, rich_est$site)

system.simpson <- t.test(filter(rich_est, reef_system=="FLA")$Simpson, filter(rich_est, reef_system=="STT")$Simpson)
system.shannon <- t.test(filter(rich_est, reef_system=="FLA")$Shannon, filter(rich_est, reef_system=="STT")$Shannon)
system.observed <- t.test(filter(rich_est, reef_system=="FLA")$Observed, filter(rich_est, reef_system=="STT")$Observed)
