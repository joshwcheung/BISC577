source("https://bioconductor.org/biocLite.R")
biocLite()
library(DNAshapeR)
library(caret)

path <- "D:/GitHub/BISC577/gcPBM/"
fa.files <- dir(path, pattern="\\.txt.fa$")
txt.files <- dir(path, pattern="\\.txt$")

## Predict DNA shapes
fn_fasta <- "D:/GitHub/BISC577/gcPBM/Mad.txt.fa"
pred <- getShape(fn_fasta)

## Encode feature vectors
##M
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType)
head(featureVector)
