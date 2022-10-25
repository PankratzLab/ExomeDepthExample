library(ExomeDepth)



file_list = list.files("~/Downloads" , pattern = ".bam$", full.names = T)

# Read in regions list for sureselect hg19 file.

exons.hg19 <-
  read.delim(
    "~/git/ExomeDepthExample/SureSelect_v7_hg19.bed.gz",
    sep = '\t',
    header = F,
    stringsAsFactors = F
  )
colnames(exons.hg19) <- c('chromosome', 'start', 'end')

# Load bam files.

my.counts <-
  getBamCounts(bed.frame = exons.hg19, bam.files = file_list)

my.counts$chromosome <-
  gsub(as.character(my.counts$chromosome),
       pattern = 'chr',
       replacement = '')

# print(table(my.counts$==0))

for (col in colnames(my.counts)) {
  print(paste0(col))
  print(table(my.counts[, col] == 0))
}

# control:        4613
