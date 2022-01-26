## Code for testing effect of substrate on microbial composition

library(phyloseq)
library(vegan)
library(broom)
library(here)
source(here("code/plotting_helpers.R"))

# Loading files
taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
# subset just VI reefs
ps.stt <- ps %>% subset_samples(reef_system=="STT" & substrate.class!="Unknown") %>% prune_taxa(taxa_sums(.)>0,.)
# change to relative abundance
ps.rel.stt <- ps.stt %>% transform_sample_counts(function(x) 100*x/sum(x))
# make new variable whether substrate was live coral or not
sample_data(ps.rel.stt)$liveCoral <- sample_data(ps.rel.stt)$substrate.class =="Live coral"
# make factors
site <- factor(as.character(sample_data(ps.rel)$site), levels = c("FLA_015", "FLA_049", "FLA_073", "Flat Cay", "Black Point", "Brewer's Bay"))
reef_system <- as.character(sample_data(ps.rel)$reef_system)
# use Bray-Curtis dissimilarity
dist.stt <- distance(ps.rel.stt, method="bray")

# adonis setup
dat <- as(sample_data(ps.rel.stt), "data.frame")
perm <- how(nperm=999)
# All samples at VI, ungrouped
perm.subst <- adonis2(dist.stt ~ sample_data(ps.rel.stt)$substrate.class, dat=dat, permutations=perm) # all substrates
perm.coral <- adonis2(dist.stt ~ sample_data(ps.rel.stt)$liveCoral, dat=dat, permutations=perm) # just live coral or not

# Site level
setBlocks(perm) <- with(dat, site)
perm.site.subst <- adonis2(dist.stt ~ sample_data(ps.rel.stt)$substrate.class, dat=dat, permutations = perm)
perm.site.coral <- adonis2(dist.stt ~ sample_data(ps.rel.stt)$liveCoral, dat=dat, permutations=perm)

# check dispersions
disper.subst <- permutest(betadisper(d = dist.stt, group=dat$substrate.class))
disper.coral <- permutest(betadisper(d = dist.stt, group=dat$liveCoral))

# Tukey HSD for dispersions
subst <- sample_data(ps.rel.stt)$substrate.class
tuk.subst <- TukeyHSD(betadisper(d=dist.stt,group=subst))

# no differences detected by anosim
anosim.coral <- anosim(distance(ps.rel.stt, method='bray'), sample_data(ps.rel.stt)$liveCoral)
anosim.subst <- anosim(distance(ps.rel.stt, method='bray'), sample_data(ps.rel.stt)$substrate.class)

# save results as res object
res <- list("coral" = list("ungrouped" = list("perm.f" = perm.coral[["F"]][1],
                                              "perm.p" = perm.coral[["Pr(>F)"]][1],
                                              "perm.r" = perm.coral[["R2"]][1],
                                              "disp.f" = disper.coral$tab$F[1],
                                              "disp.p" = disper.coral$tab$`Pr(>F)`[1]),
                           "site" = list("perm.f" = perm.site.coral[["F"]][1],
                                         "perm.p" = perm.site.coral[["Pr(>F)"]][1],
                                         "perm.r" = perm.site.coral[["R2"]][1])),
            "subst" = list("ungrouped" = list("perm.f" = perm.subst[["F"]][1],
                                              "perm.p" = perm.subst[["Pr(>F)"]][1],
                                              "perm.r" = perm.subst[["R2"]][1],
                                              "disp.f" = disper.subst$tab$F[1],
                                              "disp.p" = disper.subst$tab$`Pr(>F)`[1]),
                           "site" = list("perm.f" = perm.site.subst[["F"]][1],
                                         "perm.p" = perm.site.subst[["Pr(>F)"]][1],
                                         "perm.r" = perm.site.subst[["R2"]][1])))
# save results as .rds file
saveRDS(res, file=here("processed/disper_adonis_substrate.rds"))
