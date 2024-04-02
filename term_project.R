# setting up working enviroment
wd <- "~/Documents/school/Statistical Bioinformatics/term_project"
setwd(wd)
path <- paste(wd,"/20240328_095203/Fastq", sep = "")


# preprocessing
library(Rcpp)
library(dada2)

list.files(path)

## seperate forward & reverse reads
fnFs <- sort(list.files(path,pattern = "_R1_001.fastq.gz", full.names = TRUE)) # 001 = forward reads
fnRs <- sort(list.files(path,pattern = "_R2_001.fastq.gz", full.names = TRUE)) # 002 = reverse reads

sample.names <- sapply(strsplit(basename(fnFs), "_"), '[',1) # get sample names

# inpect read quality
plotQualityProfile(fnFs[47:48]) # trim last 10 nt (240 - 250)
plotQualityProfile(fnRs[47:48]) # trim last 20 nt (230 - 250)


