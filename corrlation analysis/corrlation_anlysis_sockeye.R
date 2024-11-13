library(tidyverse)
library(corrplot)
require(remotes)
library(Hmisc)
## sockeye
data <- read.csv("/home/zys/Documents/seafood-microbiome-analysis/corrlation analysis/abundance_data/abundance_sockeye.csv",
                 as.is = TRUE)
data_sockeye <- data %>% 
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))

microbes <- data_sockeye %>%
  select(c(Day, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
           Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
           Yersinia))
chem_prop <- data_sockeye %>%
  select(-c(Day, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
            Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
            Yersinia))

res <- rcorr(as.matrix(microbes), as.matrix(chem_prop), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes))

corplot.matrix <- res$r[1:10,-c(1:10)]

corrplot(corplot.matrix, diag = TRUE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = TRUE,
         type = "lower"
         tl.cex = 0.8,
         mar = c(1,1,1,1))

### coho
data_coho <- read.csv("/home/zys/Documents/seafood-microbiome-analysis/corrlation analysis/abundance_data/abundance_coho.csv",
                      as.is = TRUE)
data_coho <- data_coho %>%
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))
microbes_coho <- data_coho %>%
  select(c(Day, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
           Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
           Yersinia))
chem_prop_coho <- data_coho %>%
  select(-c(Day, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
           Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
           Yersinia))
res <- rcorr(as.matrix(microbes_coho), as.matrix(chem_prop_coho), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes_coho))

corplot.matrix <- res$r[1:10,-c(1:10)]

corrplot(corplot.matrix, diag = TRUE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = TRUE,
         type = "lower"
         tl.cex = 0.8,
         mar = c(1,1,1,1))

### control
data_control <- read.csv("/home/zys/Documents/seafood-microbiome-analysis/corrlation analysis/abundance_data/abundance_trial3_control.csv",
                         as.is = TRUE)
data_control <- data_control %>%
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))
microbes_control <- data_control %>%
  select(c(Day,Aminobacter,Bradyrhizobium,Cupriavidus,Hydrotalea,Mesorhizobium,Pseudaminobacter,
           Pseudomonas, Ramlibacter,Sphingopyxis,Yersinia))
chem_prop_control <- data_control %>%
  select(-c(Day,Aminobacter,Bradyrhizobium,Cupriavidus,Hydrotalea,Mesorhizobium,Pseudaminobacter,
            Pseudomonas, Ramlibacter,Sphingopyxis,Yersinia))
res <- rcorr(as.matrix(microbes_control), as.matrix(chem_prop_control), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes_control))
corplot.matrix <- res$r[1:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = TRUE,
         type = "lower"
         tl.cex = 0.8,
         mar = c(1,1,1,1))
### Nisin
data_nisin <- read.csv("/home/zys/Documents/seafood-microbiome-analysis/corrlation analysis/abundance_data/abundance_trial3_nisin.csv",
                         as.is = TRUE)
data_nisin <- data_nisin %>%
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))
microbes_nisin <- data_nisin %>%
  select(c(Day,Aminobacter,Bradyrhizobium,Cupriavidus,Hydrotalea,Mesorhizobium,Pseudaminobacter,
           Pseudomonas, Ramlibacter,Sphingopyxis,Yersinia))
chem_prop_nisin<- data_nisin%>%
  select(-c(Day,Aminobacter,Bradyrhizobium,Cupriavidus,Hydrotalea,Mesorhizobium,Pseudaminobacter,
            Pseudomonas, Ramlibacter,Sphingopyxis,Yersinia))
res <- rcorr(as.matrix(microbes_nisin), as.matrix(chem_prop_nisin), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes_nisin))
corplot.matrix <- res$r[1:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = TRUE,
         type = "lower",
         tl.cex = 0.8,
         mar = c(1,1,1,1))
### Microcin
data_microcin <- read.csv("/home/zys/Documents/seafood-microbiome-analysis/corrlation analysis/abundance_data/abundance_trial3_microcin.csv",
                          as.is = TRUE)
data_microcin <- data_microcin %>%
  filter_all(any_vars(!is.na(.))) %>%
  select_if(~ any(!is.na(.)))
microbes_microcin <- data_microcin %>%
  select(c(Day,Aminobacter,Bradyrhizobium,Cupriavidus,Hydrotalea,Mesorhizobium,Pseudaminobacter,
           Pseudomonas, Ramlibacter,Sphingopyxis,Yersinia))
chem_prop_microcin <- data_microcin %>%
  select(-c(Day,Aminobacter,Bradyrhizobium,Cupriavidus,Hydrotalea,Mesorhizobium,Pseudaminobacter,
            Pseudomonas, Ramlibacter,Sphingopyxis,Yersinia))
res <- rcorr(as.matrix(microbes_microcin), as.matrix(chem_prop_microcin), type = "spearman")
res.microbe <- rcorr(as.matrix(microbes_microcin))
corplot.matrix <- res$r[1:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE,
         tl.cex = 0.8,
         mar = c(1,1,1,1))
corrplot(res.microbe$r, diag = TRUE,
         type = "lower"
         tl.cex = 0.8,
         mar = c(1,1,1,1))

