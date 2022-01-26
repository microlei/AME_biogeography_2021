# code for permutational dispersion and permanova tests for secondary analysis
library(phyloseq)
library(vegan)
library(broom)
library(here)

source(here("../code/plotting_helpers.R"))

taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

f <- zCompositions::cmultRepl(t(asvs), method="CZM", label=0) %>% t()
ps <- phyloseq(otu_table(f, taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))

ps.clr <- ps %>% microbiome::transform("clr")
ps.dist.aitch <- distance(ps.clr, method="euclidean")

dat <- as(sample_data(ps.clr), "data.frame")

# Methods tests
methodPrimer <- sample_data(ps.clr)$methodPrimer
disp.primer <- betadisper(d=ps.dist.aitch, group=methodPrimer)
disper.primer <- permutest(disp.primer)

# environment tests
# removing samples which don't have data on collection depth, temperature, or reef type. And removing reef types with <2 samples
ps.env <- subset_samples(ps.clr, !is.na(collectionDepth) & !is.na(temperature) & !reefType %in% c("", "offshore", "patch reef", "wreck"))
ps.env.dist <- distance(ps.env, method="euclidean")
perm.env <- adonis2(ps.env.dist ~ study + reefType + collectionDepth + temperature, data=as(sample_data(ps.env),"data.frame"), by="margin")
perm.env <- tidy(perm.env)

saveRDS(perm.env, file=here("processed/perm_env.rds"))
saveRDS(disper.primer, file=here("processed/disper_primer.rds"))
