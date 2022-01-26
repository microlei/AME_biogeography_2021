## Code for testing group differences in microbiome at multiple spatial levels

library(phyloseq)
library(vegan)
library(broom)
library(here)
source(here("code/plotting_helpers.R"))

# loading files
taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)) # relative abundance
dist.rel <- phyloseq::distance(ps.rel, method="bray") # Bray-Curtis dissimilarity

# Making factors
site <- factor(as.character(sample_data(ps.rel)$site), levels = c("FLA_015", "FLA_049", "FLA_073", "Flat Cay", "Black Point", "Brewer's Bay"))
reef_system <- as.character(sample_data(ps.rel)$reef_system)
site.transect <- as.character(sample_data(ps.rel)$site.transect)

# Beta dispersions/permutation tests
disper.system <- permutest(betadisper(d = dist.rel, group=reef_system))
disper.site <- permutest(betadisper(d = dist.rel, group=site))
disper.transect <- permutest(betadisper(d = dist.rel, group=site.transect))

# Tukey HSD tests
tuk.site <- TukeyHSD(betadisper(d=dist.rel,group=site))
tuk.transect <- TukeyHSD(betadisper(d=dist.rel,group=site.transect))

# Adonis/PERMANOVA tests
dat <- as(sample_data(ps.rel),"data.frame")
perm <- how(nperm=999)
perm.system <- adonis2(dist.rel ~ reef_system, data = dat) # system level
setBlocks(perm) <- with(dat, reef_system)
perm.site <- adonis2(dist.rel ~ site, data = dat, permutations = perm) # site level
setBlocks(perm) <- with(dat, site)
perm.transect <- adonis2(dist.rel ~ site.transect, permutations=perm, data=dat) # transect level
# tidy all PERMANOVA results in a tibble
res.perm <- bind_rows(tidy(perm.transect), tidy(perm.site), tidy(perm.system))
colnames(res.perm) <- c("Scope", "Df", "Sum Sq", "R2", "Pseudo-F", "P")
res.perm$Scope[1] <- "Transect (within reef)"
res.perm$Scope[4] <- "Site (between reefs)"
res.perm$Scope[7] <- "Reef system (between reef regions)"
# tidy all permutest results in a tibble
res.disp <- bind_rows(disper.transect[["tab"]], disper.site[["tab"]], disper.system[["tab"]]) %>% select(`Sum Sq`, `Mean Sq`, F, `Pr(>F)`) %>% tibble() %>%
  bind_cols(tibble(Scope=c("Transect (within reef)", "Residual", "Site (between reefs)", "Residual", "Reef System (between reef regions)", "Residual")),.)
# make all results in a list
res = list("transect" = list("disper.f"=disper.transect[["tab"]][["F"]][1],
                            "disper.p"=disper.transect[["tab"]][["Pr(>F)"]][1],
                            "perm.f"=perm.transect[["F"]][1],
                            "perm.p"=perm.transect[["Pr(>F)"]][1]),
           "site" = list("disper.f"=disper.site[["tab"]][["F"]][1],
                         "disper.p"=disper.site[["tab"]][["Pr(>F)"]][1],
                         "perm.f"=perm.site[["F"]][1],
                         "perm.p"=perm.site[["Pr(>F)"]][1]),
           "system"=list("disper.f"=disper.system[["tab"]][["F"]][1],
                         "disper.p"=disper.system[["tab"]][["Pr(>F)"]][1],
                         "perm.F"=perm.system[["F"]][1],
                         "perm.p"=perm.system[["Pr(>F)"]][1]))
# save results as an rds and write table to csv
saveRDS(res, file=here("processed/disper_adonis_tests.rds"))
readr::write_delim(res.perm, file=here("figures/table_permanova.csv"), delim=",")
