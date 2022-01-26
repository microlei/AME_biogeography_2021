# code for figure S3, depends on plot_corncob_spatial.R
library(phyloseq)
library(here)

source(here("code/plotting_helpers.R"))

# depends on plot_corncob_spatial.R

taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x))

abunddf <- taxa_sums(ps) %>% stack() %>% mutate(rankAbundance=row_number(-values)) %>% select(-values)
prevdf <- apply(X = otu_table(ps), MARGIN=1,FUN = function(x) {sum(x>0)}) %>% stack() %>% mutate(rankPrevalence=row_number(-values)) %>% select(-values)
vardf <- apply(X=otu_table(ps.rel), MARGIN=1, FUN=function(x){var(x)}) %>% stack() %>% mutate(rankVariance=row_number(-values)) %>% select(-values)

df <- left_join(abunddf, prevdf) %>% left_join(vardf)
sig_taxa <- readRDS(here("processed/corncob_all_sites.rds"))
sig_taxa <- sig_taxa$significant_taxa
df$sig <- df$ind %in% sig_taxa

ggplot(df, aes(x=rankAbundance, y=rankVariance)) + geom_point(aes(color=sig, alpha=sig)) + scale_x_log10() + scale_y_log10() +
  labs(x="Rank Abundance", y="Rank Variance") +
  theme(legend.position="none")

ggsave(filename = here("figures/ranked_taxa.pdf"), height = 81, width = 81, units = "mm")
