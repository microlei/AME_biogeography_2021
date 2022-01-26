# plots distance decay by transect. Plot not in paper, but mantel tests are
library(ggplot2)
library(phyloseq)
library(microbiome)
library(here)

source(here("code/plotting_helpers.R"))

taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
f <- zCompositions::cmultRepl(t(asvs), method="CZM", label=0) %>% t()

ps <- phyloseq(otu_table(as.matrix(f), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))

ps.clr <- ps %>% microbiome::transform("clr")
transects <- unique(sample_data(ps.clr)$site.transect)
pairwise = data.frame()
# old code (should now use as.data.frame.table)
makePairwiseDistance <- function(x){
  d <- distance(x, method="euclidean")
  return(d %>% as.matrix() %>% as.data.frame.table(responseName="aitch_distance"))
}
for(l in transects){
  pairwise <- bind_rows(pairwise, ps.clr %>% subset_samples(site.transect==l) %>% makePairwiseDistance())
}

#appends the distance along the transect for both pairs, then takes the difference
pairwise <- pairwise %>%
  left_join(x=., y=microbiome::meta(ps.clr) %>% select(sample, site, site.transect, distance), by=c("Var1"="sample")) %>% rename(site1="site", transect1="site.transect", distVar1="distance") %>%
  left_join(x=., y=microbiome::meta(ps.clr) %>% select(sample, site, site.transect, distance), by=c("Var2"="sample")) %>% rename(site2="site", transect2="site.transect", distVar2="distance") %>%
  mutate(diff=abs(distVar1-distVar2))
pairwise <- filter(pairwise, diff!=0) #remove the ones with diff 0 as they're identical

lmresult <- summary(lm(formula = diff ~ aitch_distance, data = pairwise))
eqs <- lm_eqn(lmresult)

p <- ggplot(pairwise, aes(x=diff, y=-aitch_distance)) + geom_point() + labs(x="Geographic distance (m)", y="Similarity (-Aitchinson distance)")
p + geom_smooth(method="lm", formula = y~x, se=FALSE) +
  annotate("label", x=5,y=seq(-120,-145,length=3), label=eqs, parse=TRUE, size=4)+
  theme(axis.title = element_text(size=8), plot.title = element_text(size=8, hjust = .5))

# remake distance matrix using pairwise data
trans.aitch <- xtabs(aitch_distance ~ Var1 + Var2, data=pairwise, sparse=TRUE, addNA=TRUE) %>% as.matrix()
trans.aitch[trans.aitch==0] <- NA
geo.trans <- xtabs(diff ~ Var1 + Var2, data=pairwise, sparse=TRUE, addNA=TRUE) %>% as.matrix()
geo.trans[geo.trans==0] <- NA
# cannot do mantel test because of >50% missing values since I only had pairwise distance WITHIN not between all samples in a transect
# trans_mantel <- mantel(trans.aitch, geo.trans)

saveRDS(lmresult, file=here("processed/trans_lm.rds"))
ggsave(filename=here("figures/dist_decay_transect.png"),width=81,height=81,unit="mm")
