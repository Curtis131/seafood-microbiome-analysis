# setting up working enviroment
wd <- "/home/zys/Documents/seafood-microbiome-analysis"
setwd(wd)
path <- "/home/zys/Documents/miseq_results/20241015_184250/Fastq"


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


sample.names.r <- sapply(strsplit(basename(fnRs), "_"), '[',1)
rm_samples <- function(x,y){
  sample.names.r <- sapply(strsplit(basename(x), "_"), '[',1)
  sample.names <- sapply(strsplit(basename(y), "_"), '[',1)
  diff <- setdiff(sample.names.r,sample.names)
  if (length(diff) == 0){
    return(x)
  }
  pattern <- paste(diff,collapse = "|")
  index_to_rm <- grep(pattern,x)
  x.rm <- x[-index_to_rm]
  return(x.rm)
}

fnRs <- rm_samples(fnRs,fnFs)

if(Sys.info()['sysname'] == "Linux" || Sys.info()['sysname'] == "Darwin"){
  out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs, maxN = 0,
                       maxLen = 251, minLen = 251,
                       maxEE = c(2,2),truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = TRUE) 
} else{
  out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs, maxN = 0,
                       maxLen = 251, minLen = 251,
                       maxEE = c(2,2),truncQ = 2, rm.phix = TRUE,
                       compress = TRUE, multithread = FALSE) 
}

head(out) # should we remove samples with very few (<100) reads?

# remove samples with low reads
notvalidsamples <- which(out[,2] < 100)
sample_ids <- sample_ids[-notvalidsamples]
filtFs <- filtFs[-notvalidsamples]
filtRs <- filtRs[-notvalidsamples]
# parametric error rate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

outerrorplot.forward <- plotErrors(errF,nominalQ = TRUE)
errorplot.reverse <- plotErrors(errR,nominalQ = TRUE)

# dereplication
derepFs <- derepFastq(filtFs,verbose = TRUE)
derepRs <- derepFastq(filtRs,verbose = TRUE)

names(derepFs) <- sample.names[-notvalidsamples]
names(derepRs) <- sample.names[-notvalidsamples]

# Sample inference
## here, the samples are analyzed separately.
## we can also pool the samples before analysis, which is easier to find rare varients with low reads.
## let me know which one you think is better.
dadaFs <- dada(derepFs, err=errF, multithread = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread = TRUE)


# Merge pair read
mergers <- mergePairs(dadaFs,derepFs,dadaRs,derepRs,verbose = TRUE)

# make sequence table
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
table(nchar(getSequences(seqtab))) #make sure length of ASV do not exceed length of V4 amplicon (~464bp)
## I'm debating whether or not should I covert the ASV table to OTU table.

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread=TRUE, verbose = TRUE)
chimera.precentage <- sum(seqtab.nochim)/sum(seqtab)

# preprocessing conclusions
getN <- function(x) sum(getUniques(x))
track <- cbind(out[-notvalidsamples,], sapply(dadaFs,getN), sapply(dadaRs,getN), 
               sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR","merged", "nochim")
rownames(track) <- sample.names
write.csv(track, file = "preprocessingresult_sf2.csv")

# Assign taxonomy
## I'm using the SLIVA database
download.file(url = "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
              destfile = paste0(wd,"/silva_nr99_v138.1_train_set.fa.gz"),
              method = "curl",
              quiet = FALSE)

download.file(url = "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz",
              destfile = paste0(wd,"/silva_species_assignment_v138.1.fa.gz"),
              method = "curl",
              quiet = FALSE)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz") # we don't necessarily need to do this since I'm only interested in genus distribution.

write.csv(taxa,file = "taxa_sf2.csv")
write.csv(seqtab.nochim, file = "seqtab.nochim.csv")

# MOVING ON TO DATA ANALYSIS ##################################




