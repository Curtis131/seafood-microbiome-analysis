library(tidyverse)
library(corrplot)
require(remotes)
library(Hmisc)

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
         tl.cex = 1,
         title = "sockeye corrlations")

