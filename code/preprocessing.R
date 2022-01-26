# code for preprocessing Fl and VI output from dada2
library(here)
library(tidyverse)
library(readr)
library(data.table)

#plotting stuff
library(ggplot2)

#analysis
library(phyloseq)
library(decontam)

# reading in data
taxonomy <- read.delim(here("output/taxonomy.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
asvs <- read.delim(here("output/ASVs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
metadata <- read.delim(here("output/Spatial_metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
asvSeqs <- read.delim(here("output/ASVseqs.txt"), "\t", header=TRUE, row.names=1, check.names = FALSE)
track <- read.delim(here("output/track_reads.csv"), ",", header=TRUE, row.names = 1, check.names = FALSE)
track <- track %>% mutate(read_prop = nochim/raw) #making column prop reads retained
metadata <- metadata %>% rownames_to_column(var="row") %>% left_join(dapi, by="DAPI") %>% column_to_rownames(var="row")

# making phyloseq object
ps <- phyloseq(otu_table(as.matrix(asvs), taxa_are_rows=TRUE),
               sample_data(as.data.frame(metadata)),
               tax_table(as.matrix(taxonomy)))
ps <- merge_phyloseq(ps, sample_data(track))
ps <- subset_samples(ps, sample != "FLA_Undetermined" & sample !="STT_Undetermined")

# removing mocks
ps <- subset_samples(ps, sample != "FLA_pos-mock" & sample != "STT_mock") %>% prune_taxa(taxa_sums(.)>0,.)

# determining contaminants
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

# pruning contaminants and other stuff
ps <- prune_taxa(!contamdf.prev$contaminant, ps)
ps <- subset_taxa(ps, Kingdom!="Eukaryota")
ps <- subset_taxa(ps, Order !="Chloroplast" | is.na(Order))
ps <- subset_samples(ps, Sample_or_Control=="Sample") %>% prune_taxa(taxa_sums(.)>0,.)

# removing technical replicates
ps <- subset_samples(ps, !sample_names(ps)%in%c("FLA_45_Run2", "FLA_48_Run2", "FLA_65-N715_Run2", "FLA_54-N715_Run2", "FLA_60_Run2", "FLA_267_Run2", "FLA_279_Run2", "FLA_284-N715_Run2", "FLA_293_Run2", "FLA_299-N715_Run2", "FLA_344_Run1", "FLA_348_Run2", "FLA_351_Run2"))
ps <- subset_samples(ps, sample_sums(ps)>1000)

sample_names(ps) <- sample_data(ps)$sample #change sample names now that duplicates are gone
sample_data(ps)$transect <- factor(sample_data(ps)$transect, levels = c("1", "2", "3")) #make transect into a factor rather than number
sample_data(ps)$site <- factor(sample_data(ps)$site, levels = c("FLA_015", "FLA_049", "FLA_073", "Flat Cay", "Black Point", "Brewer's Bay"))

# Add each transect as its unique variable
sample_data(ps)$site.transect <- paste(as.character(sample_data(ps)$site),as.character(sample_data(ps)$transect)) %>%
  factor(levels = c("FLA_015 1","FLA_015 2","FLA_015 3", "FLA_049 1","FLA_049 2","FLA_049 3", "FLA_073 1","Flat Cay 1","Flat Cay 2","Flat Cay 3","Black Point 1","Black Point 2","Black Point 3","Brewer's Bay 1","Brewer's Bay 2","Brewer's Bay 3"))


# make psrelmelted
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps.rel.melted <- psmelt(ps.rel)
ps.melted <- psmelt(ps)

#re-export cleaned tsv files
new.asvs <- otu_table(ps)
new.taxonomy <- tax_table(ps)
new.metadata <- as(sample_data(ps), "data.frame")
new.asvSeqs <- asvSeqs[which(asvSeqs$asv %in% taxa_names(ps)),]

# writing the cleaned ps to rds
saveRDS(ps.rel.melted, file=here("processed/psrelmelted.rds"))
saveRDS(ps.melted, file=here("processed/psmelted.rds"))
write.table(new.asvs, file=here("processed/ASVs.txt"), sep="\t", col.names = NA, quote = FALSE)
write.table(new.taxonomy, file=here("processed/taxonomy.txt"), sep="\t", col.names = NA, quote = FALSE)
write.table(new.metadata, file=here("processed/metadata.csv"), sep=",", col.names = NA, quote=FALSE)
write.table(new.asvSeqs, file=here("processed/asvSeqs.txt"), sep="\t", col.names = NA, quote = FALSE)
