#Setup
TASKDIR = "~/Documents/UW/starita_lab/treutlein_data"
setwd("~/Documents/UW/starita_lab/treutlein_data/")
library(stats)
library(data.table)

file.copy("metadata.tsv", "metadata_copy.tsv")
btmeta <- read.table("metadata_copy.tsv")
head(btmeta)
unique(btmeta$sample)
table(btmeta$sample)
table(btmeta$batch)
# I believe that h=hour and d=day after induction.
# Batch 171122 includes _d1, _d28, and _h0
# Batch 171212 includes _d14, and _d5
# sample 180221_d35 is a polyclonal NGN-iN line from 409B2 and Sc102a1.
  # they did this test to prove that the heterogeneity is evident across cells
  # in monoclonal and polyclonal lines. 
  # maybe I exclude this data?
table(btmeta$timepoint, btmeta$sample)

