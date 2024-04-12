sample_meta <- read.csv("merge.csv", row.names = 1)
sample_meta <- sample_meta[sample_meta$treatment %in%
                             c("control","Nisin","Pediocin","Divergicin","MicrocinJ25"),]


combined_indices$time <- sample_meta$time.day.
combined_indices$treatment <- sample_meta$treatment
combined_indices$treatment <- factor(combined_indices$treatment)

# anova
aof <- function(x){
  anova(aov(x~treatment))
}
aov.result.chao1 <- apply(combined_indices$Chao1, 1, aof)
boxplot(Chao1~treatment, data = combined_indices)

anova <- aov.result %>% 
  anova(aov.result) %>%
  ggplot(aes(x = treatment, y = anova, fill = treatment))

         