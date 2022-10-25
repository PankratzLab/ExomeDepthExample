# example for calling cnvs with exome depth, and excluding related individuals for use as a reference

# install openxlsx if needed
if (!require(openxlsx)) {
  install.packages("openxlsx")
}
library(ExomeDepth)
library(openxlsx)



# bam=
# Read in regions list for sureselect hg19 file.

exons.hg19 <-
  read.csv(
    "S31285117_Regions.bed",
    sep = '\t',
    header = F,
    stringsAsFactors = F
  )
colnames(exons.hg19) <- c('chromosome', 'start', 'end')

# Load bam files.

my.counts <-
  getBamCounts(bed.frame = exons.hg19, bam.files = file_list)

# Save the dataframe.

my.counts$chromosome <-
  gsub(as.character(my.counts$chromosome),
       pattern = 'chr',
       replacement = '')
