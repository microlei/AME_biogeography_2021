# code for analyzing DAPI (not included in paper, do not upload to github)

library(ggplot2)
library(here0)

source(here("code/plotting_helpers.R"))

metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)

dapi <- metadata %>% filter(!is.na(bact.ml)) %>% select(sample, site.transect, site, bact.ml, dsRed.ml)
dapi.t <- pairwise.t.test(dapi$bact.ml, dapi$site.transect)

ggplot(dapi, aes(x=site.transect, y=bact.ml, color=factor(site, levels=c("Flat Cay", "Black Point", "Brewer's Bay")))) + geom_boxplot() +
  scale_color_manual(values=cols.site[4:6]) +
  labs(y="DAPI counts of total bacteria per mL", x="Site and transect number") +
  coord_flip()+
  scale_x_discrete(limits=rev) +
  scale_y_continuous(labels = function(x) format(x, scientific=TRUE)) +
  theme(legend.position="none", axis.text = element_text(size=7), axis.title = element_text(size=8), axis.text.x.bottom=element_text(angle=45, vjust = .5)) +
  geom_text(data=tibble(x=c(7,8,9), y=rep(3.5e+06,3), label=c("a", "ab", "b")), aes(x=x, y=y, label=rev(label)), color="black")


ggsave("figures/DAPI.pdf", width=81, height=81, units="mm")
