# example for calling cnvs with exome depth, and excluding related individuals for use as a reference

# install openxlsx if needed
if (!require(openxlsx)) {
  install.packages("openxlsx")
}
library(ExomeDepth)
library(openxlsx)

# load a counts data frame, change path as needed
load("/Users/Kitty/tmp/CountsAll.RData")

# we limit this analysis to the autosome. To call chrX/chrY, simply select a single sex reference set that matches the sample of interest.
all.counts = all.counts[which(!all.counts$chromosome %in% c("chrX", "chrY")), ]

# load the family codes xlsx provided
familyCodes = read.xlsx("/Users/Kitty/Downloads/Family Codes for CNV .xlsx")

# Add and "X" to the id if the first character of the id is a number.
# This allows us to match columns of all.counts to the Bam.ID with a new column named Bam.ID
familyCodes$all.counts.ID = familyCodes$Bam.ID
# These IDs need to be updated since they start with an number
update = grepl("^[[:digit:]]+", familyCodes$all.counts.ID)
#Update them
familyCodes$all.counts.ID[update] = paste0("X", familyCodes$all.counts.ID[update])

# all.cnvs will store all of our CNV calls
all.cnvs = data.frame()
for (id in familyCodes$all.counts.ID) {
  if (!id %in% colnames(all.counts)) {
    print(paste0("WARNING: ", id, " was not found in the all.counts dataframe"))
    next
  }
  
  # The unique family ID
  currentFamily = unique(familyCodes$Family.ID[which(familyCodes$all.counts.ID ==
                                                       id)])
  # Exclude all members of the family from the reference set
  excludeFromReferenceSet = familyCodes$all.counts.ID[which(familyCodes$Family.ID ==
                                                              currentFamily)]
  
  # print a message regardign which samples are excluded
  print(
    paste0(
      "For id ",
      id,
      " the following ",
      length(excludeFromReferenceSet),
      " members of family ",
      currentFamily,
      " will be excluded from the reference set (including self): ",
      paste0(excludeFromReferenceSet, collapse = ", ")
    )
  )
  # This is the sample we are calling CNVS for
  my.test <- all.counts[, id]
  
  # All possible samples, i.e all samples in the counts data frame...i.e all columns that are not chrom
  allPossibleReferenceSamples <-
    colnames(all.counts)[!colnames(all.counts) %in% c("chromosome", "start", "end", "exon")]
  
  # in case there are duplicated sample names, get a unique list of samples
  allPossibleReferenceSamples = unique(allPossibleReferenceSamples)
  
  # Reference samples (those samples not in the same family). i.e remove samples in the excludeFromReferenceSet
  my.ref.samples <-
    allPossibleReferenceSamples[!allPossibleReferenceSamples %in% excludeFromReferenceSet]
  
  print(paste0(
    "available reference samples goes from ",
    length(allPossibleReferenceSamples),
    " to ",
    length(my.ref.samples)
  ))
  if (length(allPossibleReferenceSamples) - length(my.ref.samples) != length(excludeFromReferenceSet)) {
    miss = excludeFromReferenceSet[!excludeFromReferenceSet %in% allPossibleReferenceSamples]
    print(
      paste0(
        "warning, assuming that some samples from family ",
        currentFamily,
        " were not available in the counts data frame - these ",
        length(miss),
        " were missing for exclusion ",
        paste0(miss, collapse = ", ")
      )
    )
  }
  
  # we now extract the counts for all reference samples that are not related to our sample of interest
  my.reference.set <- as.matrix(all.counts[, my.ref.samples])
  
  # from the unrelated reference samples, select an optimized set
  my.choice <- select.reference.set (
    test.counts = my.test,
    reference.counts = my.reference.set,
    bin.length = (all.counts$end - all.counts$start) /
      10,
    n.bins.reduced = 10000
  )
  
  # extract the read counts for the optimzed set of reference samples
  my.matrix <-
    as.matrix(all.counts[, my.choice$reference.choice, drop = FALSE])
  
  
  # for each exon, sum the counts across each row
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum)
  all.exons <- new(
    'ExomeDepth',
    test = my.test,
    reference = my.reference.selected,
    formula = 'cbind(test, reference) ~ 1'
  )
  
  # call cnvs
  all.exons <- CallCNVs(
    x = all.exons,
    transition.probability = 10 ^ -4,
    chromosome = all.counts$chromosome,
    start = all.counts$start,
    end = all.counts$end,
    name = all.counts$exon
  )
  
  CNV_calls <- all.exons@CNV.calls
  
  # The following additions to dataframe can be changed as desired
  
  # store the sample id in this dataframe
  CNV_calls$sample = id
  CNV_calls$familyID = familyID
  # store the samples that were selected for the reference distribution
  CNV_calls$referenceSelected = paste0(my.choice$reference.choice, collapse = ",")
  
  # store the samples that were not considered for the reference distribution
  CNV_calls$referenceExcluded = paste0(excludeFromReferenceSet, collapse = ",")
  
  # combine this current samples calls with all others
  all.cnvs = rbind(all.cnvs, CNV_calls)
  
}

# write the calls to a file
write.table(
  all.cnvs,
  file = "all.cnvs.txt",
  quote = F,
  col.names = T,
  sep = "\t",
  row.names = F
)

#write all.cnvs to file if desired