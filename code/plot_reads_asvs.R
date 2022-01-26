# code for figure S4
library(phyloseq)
library(here)

source(here("code/plotting_helpers.R"))

# current study
taxonomy <- read.delim(here("processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
sample_data(ps)$study <- "Present study"

# meta-analysis
taxonomy <- read.delim(here("SMeta/processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("SMeta/processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("SMeta/processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

ps.meta <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows = TRUE),
                    sample_data(as.data.frame(metadata)),
                    tax_table(as.matrix(taxonomy)))

dt <- tibble(study = c(sample_data(ps.meta)$study, sample_data(ps)$study),
                 reads = c(colSums(otu_table(ps.meta)), colSums(otu_table(ps))),
                 ASVs = c(colSums(otu_table(ps.meta)!=0), colSums(otu_table(ps)!=0)))

p <- ggplot(dt, aes(x=reads, y=ASVs)) + geom_point(aes(color=study)) +
  scale_color_manual(name="Study", labels=c(labels.study, "Present study"), values=c(rep("grey",5),"blue")) +
  scale_x_log10() +
  labs(x="Total reads", y="Observed ASVs") +
  theme(axis.title = element_text(size=8), legend.position="none")

p + geom_smooth(data=filter(dt, study!="Present study"), method=lm, se=FALSE, color="black") +
  geom_smooth(data=filter(dt, study=="Present study"), method=lm, se=FALSE, color="red") +

ggsave(file=here("figures/reads_asvs.pdf"), width=81, height=81, units = "mm")
