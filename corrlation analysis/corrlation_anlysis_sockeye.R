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
           Vibrio, Yersinia))
chem_prop <- data_sockeye %>%
  select(-c(Day, Brochothrix, Carnobacterium, Clostridium.sensu.stricto.1,
            Cupriavidus, Latilactobacillus, Leuconostoc, Photobacterium, Pseudomonas,
            Vibrio, Yersinia))
colnames(chem_prop) <- gsub("^X|Average_", "", colnames(chem_prop))


res <- rcorr(as.matrix(microbes), as.matrix(chem_prop), type = "pearson")
res.microbe <- rcorr(as.matrix(microbes))
cor_test <-res$P[2:11, -c(1:11)]
cor_microbe_test <- res.microbe$P

corplot.matrix <- res$r[2:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE, method = "square",
         p.mat = cor_test, 
         tl.col = "black",
         addCoef.col = "black",
         cl.cex = 0.4,
         number.cex = 0.4,
         tl.cex = 0.6, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "yellow",
         pch.cex = 1)

corrplot(res.microbe$r, diag = TRUE, method = "square",
         type = "lower",
         tl.col = "black",
         p.mat = cor_microbe_test,
         addCoef.col = "black",
         number.cex = 0.4,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.col = "yellow",
         pch.cex = 1,
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
           Vibrio, Yersinia))
colnames(chem_prop_coho) <- gsub("^X|Average_", "", colnames(chem_prop_coho))


res <- rcorr(as.matrix(microbes_coho), as.matrix(chem_prop_coho), type = "pearson")
res.microbe <- rcorr(as.matrix(microbes_coho))
cor_test <-res$P[2:11, -c(1:11)]
cor_microbe_test <- res.microbe$P

corplot.matrix <- res$r[2:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE, method = "square",
         p.mat = cor_test, 
         tl.col = "black",
         addCoef.col = "black",
         cl.cex = 0.4,
         number.cex = 0.4,
         tl.cex = 0.6, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "yellow",
         pch.cex = 1)

corrplot(res.microbe$r, diag = TRUE, method = "square",
         type = "lower",
         tl.col = "black",
         p.mat = cor_microbe_test,
         addCoef.col = "black",
         number.cex = 0.4,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.col = "yellow",
         pch.cex = 1,
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
colnames(chem_prop_control) <- gsub("^X|Average_", "", colnames(chem_prop_control))


res <- rcorr(as.matrix(microbes_control), as.matrix(chem_prop_control), type = "pearson")
res.microbe <- rcorr(as.matrix(microbes_control))
cor_test <-res$P[2:11, -c(1:11)]
cor_microbe_test <- res.microbe$P

corplot.matrix <- res$r[2:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE, method = "square",
         p.mat = cor_test, 
         tl.col = "black",
         addCoef.col = "black",
         cl.cex = 0.4,
         number.cex = 0.4,
         tl.cex = 0.6, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "yellow",
         pch.cex = 1)

corrplot(res.microbe$r, diag = TRUE, method = "square",
         type = "lower",
         tl.col = "black",
         p.mat = cor_microbe_test,
         addCoef.col = "black",
         number.cex = 0.4,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.col = "yellow",
         pch.cex = 1,
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
colnames(chem_prop_nisin) <- gsub("^X|Average_", "", colnames(chem_prop_nisin))


res <- rcorr(as.matrix(microbes_nisin), as.matrix(chem_prop_nisin), type = "pearson")
res.microbe <- rcorr(as.matrix(microbes_nisin))
cor_test <-res$P[2:11, -c(1:11)]
cor_microbe_test <- res.microbe$P

corplot.matrix <- res$r[2:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE, method = "square",
         p.mat = cor_test, 
         tl.col = "black",
         addCoef.col = "black",
         cl.cex = 0.4,
         number.cex = 0.4,
         tl.cex = 0.6, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "yellow",
         pch.cex = 1)

corrplot(res.microbe$r, diag = TRUE, method = "square",
         type = "lower",
         tl.col = "black",
         p.mat = cor_microbe_test,
         addCoef.col = "black",
         number.cex = 0.4,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.col = "yellow",
         pch.cex = 1,
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
colnames(chem_prop_microcin) <- gsub("^X|Average_", "", colnames(chem_prop_microcin))


res <- rcorr(as.matrix(microbes_microcin), as.matrix(chem_prop_microcin), type = "pearson")
res.microbe <- rcorr(as.matrix(microbes_microcin))
cor_test <-res$P[2:11, -c(1:11)]
cor_microbe_test <- res.microbe$P

corplot.matrix <- res$r[2:11,-c(1:11)]

corrplot(corplot.matrix, diag = TRUE, method = "square",
         p.mat = cor_test, 
         tl.col = "black",
         addCoef.col = "black",
         cl.cex = 0.4,
         number.cex = 0.4,
         tl.cex = 0.6, 
         insig = "label_sig",
         sig.level = c(0.001, 0.01, 0.05),
         pch.col = "yellow",
         pch.cex = 1)

corrplot(res.microbe$r, diag = TRUE, method = "square",
         type = "lower",
         tl.col = "black",
         p.mat = cor_microbe_test,
         addCoef.col = "black",
         number.cex = 0.4,
         sig.level = c(0.001, 0.01, 0.05),
         insig = "label_sig",
         pch.col = "yellow",
         pch.cex = 1,
         tl.cex = 0.8,
         mar = c(1,1,1,1))

