# code for figure S1
source(here("code/plotting_helpers.R"))
library(cowplot)

metadata <- read.delim(here("processed/metadata.csv"), ",", header = TRUE, row.names=1, check.names=FALSE)
benthos_vi <- metadata %>% filter(reef_system=="STT") %>% select(sample, site, transect, distance, substrate.class, site.transect)
benthos_vi <- benthos_vi %>% group_by(site) %>% count(substrate.class) %>% mutate(cover=100*n/sum(n)) %>% group_by(site)
benthos_fl <- read.delim(here("processed/fl_benthos.csv"), ",", header=TRUE)
colnames(benthos_fl) <- c("site", "substrate.class", "cover")

p_fl <- ggplot(benthos_fl, aes(x=site, y=cover, fill=substrate.class)) + geom_bar(stat="identity", position="stack") + labs(fill="Substrate class", x="", y="% cover") +
  theme(axis.text.x.bottom = element_text(angle=40, vjust=.7))
p_vi <- ggplot(benthos_vi, aes(x=site, y=cover, fill=substrate.class)) + geom_bar(stat="identity", position="stack") + labs(fill="Substrate class", x="", y="% cover") +
  theme(axis.text.x.bottom = element_text(angle=40, vjust =.7))

plot_grid(p_fl, p_vi, labels="AUTO")
ggsave(filename=here("figures/substrate.pdf"), width=160, height=100, units="mm")
readr::write_delim(benthos_vi, file=here("processed/vi_benthos.csv"), delim=",")
