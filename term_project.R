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

# filter & trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## there are 4 samples that only have reverse reads but no forward reads.
## I'll check if the data didn't upload properly or there weren't any forward reads at the first place.
## For now, I'll remove the 4 samples.
sample.names.r <- sapply(strsplit(basename(fnRs), "_"), '[',1)
rm_samples <- function(x,y){
  sample.names.r <- sapply(strsplit(basename(x), "_"), '[',1)
  sample.names <- sapply(strsplit(basename(y), "_"), '[',1)
  diff <- setdiff(sample.names.r,sample.names)
  pattern <- paste(diff,collapse = "|")
  index_to_rm <- grep(pattern,x)
  x.rm <- x[-index_to_rm]
  return(x.rm)
}

fnRs <- rm_samples(fnRs,fnFs)

if(Sys.info()['sysname'] == "Linux" || Sys.info()['sysname'] == "macOs"){
  out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs, truncLen = c(240,230), maxN = 0,
                       maxEE = c(2,2),truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = TRUE) 
} else{
  out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs, truncLen = c(240,230), maxN = 0,
                       maxEE = c(2,2),truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = FALSE) 
}

head(out)
