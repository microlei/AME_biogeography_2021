# Code for corncob differential abundance tests. Figures not used in paper, but differentially abundant taxa are recorded here
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(corncob)
library(here)
source("code/plotting_helpers.R")

# loading files
taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
asvSeqs <- read.delim(here("processed/asvSeqs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))

# Removing any taxa with fewer than 10 reads total
ps.abund <- prune_taxa(taxa_sums(ps)>10,ps)
ps.abund <- add_taxonomy_column(ps.abund)
sample_data(ps.abund)$site <- factor(sample_data(ps.abund)$site, levels=breaks.site)

da.site <- differentialTest(formula = ~ site,
                            phi.formula = ~ site,
                            formula_null = ~ 1,
                            phi.formula_null = ~ site,
                            test="Wald", boot=FALSE,
                            data=ps.abund,
                            fdr_cutoff=0.05)
da.data <- cleanDA(da.site)
da.data <- da.data %>% mutate(tax=paste0(Taxonomy," (", ASV,")"))

ggplot(da.data %>% filter(site=="Brewer's Bay" & p<0.05), aes(x=Estimate, y=sort(tax, decreasing=FALSE))) +
  geom_vline(xintercept = 0, color="gray50", lty="dashed", alpha=0.75, lwd=1)+
  geom_point()+
  geom_errorbarh(aes(xmin=Estimate-1.96*StdE, xmax=Estimate+1.96*StdE), height=.3)+
  theme_light() + labs(title="Brewer's Bay")

ggplot(da.data %>% filter(site=="Flat Cay" & p<0.05), aes(x=Estimate, y=sort(tax, decreasing=FALSE))) +
  geom_vline(xintercept = 0, color="gray50", lty="dashed", alpha=0.75, lwd=1)+
  geom_point()+
  geom_errorbarh(aes(xmin=Estimate-1.96*StdE, xmax=Estimate+1.96*StdE), height=.3)+
  theme_light() + labs(title="Flat Cay")

ggplot(da.data %>% filter(site=="FLA_049" & p<0.05), aes(x=Estimate, y=sort(tax, decreasing=FALSE))) +
  geom_vline(xintercept = 0, color="gray50", lty="dashed", alpha=0.75, lwd=1)+
  geom_point()+
  geom_errorbarh(aes(xmin=Estimate-1.96*StdE, xmax=Estimate+1.96*StdE), height=.3)+
  theme_light() + labs(title="Site 49")

ggplot(da.data %>% filter(site=="FLA_073" & p<0.05), aes(x=Estimate, y=sort(tax, decreasing=FALSE))) +
  geom_vline(xintercept = 0, color="gray50", lty="dashed", alpha=0.75, lwd=1)+
  geom_point()+
  geom_errorbarh(aes(xmin=Estimate-1.96*StdE, xmax=Estimate+1.96*StdE), height=.3)+
  theme_light() + labs(title="Site 73")

p <- ggplot(da.data, aes(x=x, y=sort(taxa, decreasing = TRUE))) +
  geom_vline(xintercept = 0, color = "gray50", lty = "dashed", alpha = 0.75, lwd = 1) +
  geom_point() +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = .3) +
  theme_light() +
  facet_wrap(~variable, scales = "free_x", nrow = 1, labeller = labeller(variable=site.labs)) +
  labs(title = "", x = "", y = "Taxa") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(da.site, here("processed/corncob_all_sites.rds"))
sig.taxa <- bind_cols(taxonomy, asvSeqs) %>% select(-asv)
sig.taxa <- sig.taxa[da.site$significant_taxa,]
write.csv(sig.taxa, file=here("figures/corncob_spatial_sig_taxa.csv"))
ggsave(p, filename = here("figures/corncob.pdf"))
