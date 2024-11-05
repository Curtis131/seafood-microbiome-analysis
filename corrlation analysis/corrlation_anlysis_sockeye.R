library(tidyverse)
library(corrplot)
require(remotes)
library(Hmisc)
## sockeye
data <- read.csv("/home/zys/Downloads/abundance_data/abundance_sockeye_corr.csv",
                 as.is = TRUE)
data_sockeye <- data %>% 
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))

microbes <- data_sockeye %>%
  select(c(Day,Bradyrhizobium, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
           Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
           Yersinia))
chem_prop <- data_sockeye %>%
  select(-c(Day,Bradyrhizobium, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
            Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
            Yersinia))

res <- rcorr(as.matrix(microbes), as.matrix(chem_prop), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes))

corplot.matrix <- res$r[1:11,-c(1:11)]

corrplot(corplot.matrix, diag = FALSE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = FALSE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))

### coho
data_coho <- read.csv("/home/zys/Documents/seafood-microbiome-analysis/corrlation analysis/abundance_data/abundance_corr_coho.csv",
                      as.is = TRUE)
data_coho <- data_coho %>%
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))
microbes_coho <- data_coho %>%
  select(c(Day,Bradyrhizobium, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
           Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
           Yersinia))
chem_prop_coho <- data_coho %>%
  select(-c(Day,Bradyrhizobium, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
           Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
           Yersinia))
res <- rcorr(as.matrix(microbes_coho), as.matrix(chem_prop_coho), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes_coho))

corplot.matrix <- res$r[1:11,-c(1:11)]

corrplot(corplot.matrix, diag = FALSE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = FALSE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))

### 