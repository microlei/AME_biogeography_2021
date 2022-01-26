# code for testing effect of substrate on Fl and VI composition
library(here)
library(phyloseq)
library(tidyverse)
library(corncob)

taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
ps.stt <- ps %>% subset_samples(reef_system=="STT" & substrate.class!="Unknown") %>% prune_taxa(taxa_sums(.)>0,.)
sample_data(ps.stt)$liveCoral <- sample_data(ps.stt)$substrate.class =="Live coral"

# Corncob detects no differentially abundant taxa between Live Coral and other substrate types
# abundant taxa used to save computational time, but results are the same with the full otu table
set.seed(1)
ps.abund.stt <- prune_taxa(taxa_names(ps.rel.stt %>% prune_taxa(taxa_sums(.)/nsamples(.) > 0.001,.)),ps)
bbdml(formula= ASV1 ~ site, phi.formula = ~ site, data=ps.abund.stt)
da.coral.abund <- differentialTest(formula = ~ liveCoral,
                                   phi.formula = ~ liveCoral,
                                   formula_null = ~ 1,
                                   phi.formula_null = ~ liveCoral,
                                   test="Wald", boot=FALSE,
                                   data=ps.abund.stt,
                                   fdr_cutoff=0.05)
da.coral.abund$significant_taxa
