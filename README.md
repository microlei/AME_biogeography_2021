# Biogeography of reef water microbes from within reef to global scales

Code used for analyses, figure generation, and formatting of the manuscript.

## Sequence quality control

Raw sequence reads for the primary analysis were processed using the DADA2 workflow optimized for iSeq data found in [this](https://github.com/microlei/apprill-iseq) repository, using parameters detailed in the paper. Output from DADA2 are found in the folder 'output'. For the secondary analysis, the MiSeq version found [here](https://github.com/microlei/apprill-miseq) was used and output can be found in 'SMeta/output'.

Reads were then filtered and preprocessed using 'preprocessing.R' and 'preprocessing_SMeta.R' and deposited in the respective 'processed' folder.

## Data analysis

All code for analyses, figures, and tables generated for this paper (and some that didn't make it) can be found in the 'code' folder with descriptive filenames and a short descriptive header indicating the number of the figure it generates.

## Manuscript

The manuscript for this paper was written in RMarkdown in the file 'Spatial_manuscript.Rmd' and exported into Word format. To compile the manuscript after running all scripts in 'code', run the following in R:

```R
library(rmarkdown)
library(bookdown)
library(tools)
render("Spatial_manuscript.Rmd", output_format="word_document2")
```
