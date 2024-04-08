# Script used to generate experiment metadata 1 $ 2

## csv bactericocin treatments
wd <- "~/Documents/school/Statistical Bioinformatics/term_project"
setwd(wd)
path <- paste(wd,"/20240328_095203/Fastq", sep = "")

fnFs <- sort(list.files(path,pattern = "_R1_001.fastq.gz", full.names = TRUE))
sample.name <- sapply(strsplit(basename(fnFs), "_"), '[',1)
sample.name <- sort(sample.name[125:198])

time <- sapply(strsplit(sample.name, "-"), '[',2)
treatment <- sapply(strsplit(sample.name, "-"), '[',1)
treatment <- sapply(strsplit(sample.name, ""), '[',1)
treatment <- factor(treatment,levels = c("A","B","C","D","E"),
       labels = c("control","Nisin","Pediocin","Divergicin","MicrocinJ25"))


bacterocin.treatment <- data.frame(sample.name,time,treatment)
