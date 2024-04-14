library(ggplot2)
library(dplyr)


combined_indices <- read.csv("combined_diversity_indices.csv", row.names = 1)
sample_meta <- read.csv("merge.csv", row.names = 1)
sample_meta <- sample_meta[sample_meta$treatment %in%
                             c("control","Nisin","Pediocin","Divergicin","MicrocinJ25"),]



combined_indices$time <- bacterocin_exp@sam_data$time.day.
combined_indices$treatment <- bacterocin_exp@sam_data$treatment
combined_indices$treatment <- factor(combined_indices$treatment)
combined_indices$time <- factor(combined_indices$time)

# take average of replicates
combined_indices <- combined_indices %>%
  group_by(time,treatment) %>%
  summarise(
    Chao1 = mean(Chao1),
    Shannon = mean(Shannon),
    Simpson = mean(Simpson),
    Observed = mean(Observed)
  )

# log transformation
combined_indices <- combined_indices %>%
  mutate( Chao1 = log(Chao1),
          Shannon = log(Shannon) + 0.1590,
          Simpson = log(Simpson) + 0.158972,
          Observed = log(Observed))


# anova chao1
aof <- function(x){
  y <- anova(aov(x ~ treatment, data = combined_indices))
  return(y)
}

aof.plot <- function(x,y){
  plot <- ggplot(combined_indices, aes(x = treatment, y = x, fill = treatment)) +
  geom_boxplot() +
    # facet_wrap( ~ time, ncol = 2) +
    theme_classic() + 
    ylab(colnames(combined_indices[y])) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(plot)
}

anova_chao1 <- aof(combined_indices$Chao1)

aof.plot(combined_indices$Chao1,3)

# anova shannon
anova_shannon <- aof(combined_indices$Shannon)
aof.plot(combined_indices$Shannon,4)

# anova simpson
anova_simpson <- aof(combined_indices$Simpson)
aof.plot(combined_indices$Simpson,5)

# anova 
anova_observed <- aof(combined_indices$Observed)
aof.plot(combined_indices$Observed, 6)

