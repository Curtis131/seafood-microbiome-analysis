library(tibble)
library(dplyr)

sample.dat <- read.csv("preprocessingresult_sf2.csv",row.name = 1)
metadat <- as_tibble(rownames(sample.dat))
metadat <- metadat %>%
  rename(sample = value) %>%
  add_column(treatment = 0, day = 0, replicate = 0) %>%
  mutate(
    treatment = sapply(strsplit(sample, "-"), '[', 1),
    day = sapply(strsplit(sample, "-"), '[', 2),
    # Extracting the last two characters for replicate
    replicate = substr(sample, nchar(sample) - 2, nchar(sample))
  )
         
write.csv(metadat, file = "metadatasf2.csv")
