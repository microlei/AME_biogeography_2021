# preprocessing dada2 output for secondary analysis

library(tidyverse)
library(readr)
library(data.table)
library(here)
#plotting stuff
library(ggplot2)
#analysis
library(phyloseq)

taxonomy <- read.delim(here("SMeta/output/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("SMeta/output/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("SMeta/output/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
asvSeqs <- read.delim(here("SMeta/output/asvSeqs.txt"), "\t", header=TRUE, row.names=1, check.names=FALSE)
track <- read.delim(here("SMeta/output/track.tsv"), "\t", header=TRUE, row.names=1)
track <- mutate(track, prop.retained = nochim/reads.in)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))

ps <- merge_phyloseq(ps, sample_data(track))

ps <- subset_samples(ps, nochim >10000)
ps <- subset_taxa(ps, Kingdom!="Eukaryota")
ps <- subset_taxa(ps, Order !="Chloroplast" | is.na(Order))

new.asvs <- otu_table(ps)
new.taxonomy <- tax_table(ps)
new.metadata <- as(sample_data(ps), "data.frame")
new.asvSeqs <- asvSeqs[which(asvSeqs$asv %in% taxa_names(ps)),]

write.table(new.asvs, file=(here"SMeta/processed/ASVs.txt"), sep="\t", col.names = NA, quote = FALSE)
write.table(new.taxonomy, file=here("SMeta/processed/taxonomy.txt"), sep="\t", col.names = NA, quote = FALSE)
write.table(new.metadata, file=here("SMeta/processed/metadata.csv"), sep=",", col.names = NA, quote=FALSE)
write.table(new.asvSeqs, file=here("SMeta/processed/asvSeqs.txt"), sep="\t", col.names = NA, quote = FALSE)
