# code for venn diagram overlap of secondary analysis (not included in paper)
library(VennDiagram)
library(ggplot2)
library(phyloseq)
library(here)

source(here("code/plotting_helpers.R"))

taxonomy <- read.delim(here("SMeta/processed/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("SMeta/processed/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("SMeta/processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
asvSeqs <- read.delim(here("SMeta/processed/asvSeqs.txt"), "\t", header=TRUE, row.names=1, check.names=FALSE)

ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))

taxa_list <- list(Apprill_2021 = subset_samples(ps, study=="Apprill et al 2021") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  Becker_2020 = subset_samples(ps, study=="Becker et al 2020") %>% prune_taxa(taxa_sums(.)>0,.) %>%  taxa_names(),
                  Becker_2021 = subset_samples(ps,study=="Becker et al 2021") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  Neave_2017 = subset_samples(ps, study=="Neave et al 2017") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names(),
                  Weber_2020 = subset_samples(ps, study=="Weber et al 2020") %>% prune_taxa(taxa_sums(.)>0,.) %>% taxa_names())

overlap <- calculate.overlap(taxa_list)
overlap_asvs <- bind_cols(taxonomy, asvSeqs) %>% select(-asv)
overlap_asvs <- overlap_asvs[overlap$a31,]

write.csv(overlap_asvs, file=here("figures/SMeta_overlap.csv"))
saveRDS(taxa_list, here("SMeta/processed/taxa_list.rds"))
saveRDS(overlap, file=here("SMeta/processed/SMeta_overlap.rds"))

venn.plot <- venn.diagram(taxa_list, category.names = labels.study,
                          filename = here("figures/meta_venn.png"),
                          main="Venn diagram of taxa shared between studies",
                          fontfamily="sans",
                          imagetype="png",
                          resolution=300,
                          width=105,
                          height=105,
                          units="mm",
                          cex=.75,
                          cat.cex = 0.75,
                          cat.pos=c(0,-25,-130, 180,0),
                          cat.dist=c(.2,.22,.25,.2,.2),
                          fill=cols.study[1:5])
