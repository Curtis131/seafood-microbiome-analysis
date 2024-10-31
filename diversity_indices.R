# Load necessary libraries
library(phyloseq)
library(tidyverse)
library(SRS)

# Load the ASV and taxonomy data
asv_data <- read.csv("seqtab.nochim.csv",row.names = 1)

# seperate 2 experiments for salmon & shrimp
taxonomy_data <- read.csv("taxa_sf2.csv",row.names = 1) 
metadata <- read.csv("metadatasf2.csv", row.names = 1)
salmonexp <- metadata %>%
  filter(treatment %in% c("A","B","E"))
shrimpexp <- metadata %>%
  filter(treatment %in% c("C","I","RB"))
salmon_asv <- asv_data %>%
  filter(rownames(asv_data) %in% salmonexp$sample)
shrimp_asv <- asv_data %>%
  filter(rownames(asv_data) %in% shrimpexp$sample)

# normalization
salmon_asv_srs <- salmon_asv %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin = 5000) %>%
  t() %>%
  as.matrix()
colnames(salmon_asv_srs) <- colnames(salmon_asv)

shrimp_asv_srs <- shrimp_asv %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin = 5330) %>%
  t() %>%
  as.matrix()
colnames(shrimp_asv_srs) <- colnames(shrimp_asv)

# Convert taxonomy data to matrix and specify column names for Phyloseq
taxonomy_matrix <- as.matrix(taxonomy_data)
colnames(taxonomy_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# create phyloseq objects
## shrimp ##
shrimp_otu_table <- otu_table(as.matrix(shrimp_asv_srs), taxa_are_rows = FALSE)
tax.table <- tax_table(taxonomy_matrix)
shrimpexp <- shrimpexp %>%
  filter(sample %in% rownames(shrimp_asv_srs))
  
shrimp_sample <- sample_data(shrimpexp)
rownames(shrimp_sample) <- shrimp_sample$sample

shrimp <- phyloseq(shrimp_otu_table,tax.table,shrimp_sample)

## salmon ##
salmon_otu_table <- otu_table(as.matrix(salmon_asv_srs), taxa_are_rows = FALSE)
salmonexp <- salmonexp %>%
  filter(sample %in% rownames(salmon_asv_srs))

salmon_sample <- sample_data(salmonexp)
rownames(salmon_sample) <- salmon_sample$sample

salmon <- phyloseq(salmon_otu_table,tax.table,salmon_sample)

# analysis
## shrimp ##
shrimp <- prune_taxa(taxa_sums(shrimp) >0, shrimp) # ASV with no counts
shrimp.genus <- tax_glom(shrimp, taxrank = "Genus") # agglomerate genus 
shrimp.genus.df <- psmelt(shrimp.genus)
shrimp.genus.df <- shrimp.genus.df %>%
  mutate(relative_abundance = sapply(Abundance, function(x) x/sum(Abundance)))
shrimp.genus.df$Genus <- as.factor(shrimp.genus.df$Genus)
shrimp.genus.df.ag <- shrimp.genus.df %>%
  group_by(Sample, Genus,day,treatment) %>%
  summarize(total_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(total_sample_abundance = sum(total_abundance)) %>%
  ungroup() %>%
  mutate(relative_abundance = total_abundance / total_sample_abundance) %>%
  select(Sample, day, treatment, Genus,total_abundance, relative_abundance) %>%
  filter(relative_abundance > 0) %>%
  mutate(Genus = ifelse(relative_abundance < 0.01, "others", as.character(Genus))) %>%
  group_by(Sample, Genus,treatment,day) %>%
  summarize(total_abundance = sum(total_abundance), 
            relative_abundance = sum(relative_abundance), 
            .groups = 'drop') %>%
  mutate(treatment = factor(treatment))
shrimp.genus.df.ag$day <- factor(shrimp.genus.df.ag$day,
                                   levels = c(0,2,4,6,8,10))
shrimp.genus.df.ag <- arrange(shrimp.genus.df.ag,shrimp.genus.df.ag$day)
shrimp.genus.df.ag$Sample <- factor(shrimp.genus.df.ag$Sample, 
                                          levels = unique(shrimp.genus.df.ag$Sample[order(shrimp.genus.df.ag$day)]))

######## abundance  plot ###########
ggplot(data = shrimp.genus.df.ag,
       aes(x = Sample, y = relative_abundance, fill = Genus))+
  geom_bar(aes(),stat = "identity", position = "stack") +
  facet_grid(day ~ ., scales = "free_y") +
  coord_flip()
  
# or #
shrimp.genus.df.ag <- shrimp.genus.df.ag %>%
  mutate(day_treatment = paste(day, treatment, sep = " - "))  # Combine day and treatment
  
ggplot(data = shrimp.genus.df.ag %>% subset(day == "0"),
        aes(x = Sample, y = relative_abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  facet_wrap(~ day_treatment, scales = "free_y") +  # Facet by combined variable
  labs(x = "Sample", y = "Relative Abundance", fill = "Genus") +
  theme_minimal()  # Optional: use a minimal theme for better clarity

################
  
# Calculate Chao1 diversity
chao1_shrimp <- estimate_richness(shrimp, measures = "Chao1")
shannon_shrimp <- estimate_richness(shrimp, measures = "Shannon")
simpson_shrimp <- estimate_richness(shrimp, measures = "Simpson")
observed_shrimp <- estimate_richness(shrimp, measures = "Observed")

combined_indices_shrimp <- data.frame(Sample = rownames(shannon_shrimp),
                               Chao1 = chao1_shrimp$Chao1,
                               Shannon = shannon_shrimp$Shannon,
                               Simpson = simpson_shrimp$Simpson, 
                               Observed = observed_shrimp$Observed,
                               day = shrimp@sam_data$day,
                               replicate = shrimp@sam_data$replicate,
                               treatment = shrimp@sam_data$treatment)

shrimp@sam_data$day <- factor(shrimp@sam_data$day,
                              levels = c(0,2,4,6,8,10))
plot_richness(shrimp, x = "day", measures = "Shannon") +
  geom_boxplot() +
  ylab("Shannon")

plot_richness(shrimp, x = "day", measures = "Observed") +
  geom_boxplot() +
  ylab("Observed")

plot_richness(shrimp, x = "day", measures = "Simpson") +
  ylab("Simpson") +
  geom_boxplot()

plot_richness(shrimp, x = "day", measures = "Chao1") +
  ylab("Chao1") +
  geom_boxplot()

# take average of replicates, log transformation
combined_indices_shrimp_trans <- combined_indices_shrimp %>%
  group_by(day,treatment) %>%
  summarise(
    Chao1 = mean(Chao1),
    Shannon = mean(Shannon),
    Simpson = mean(Simpson),
    Observed = mean(Observed)
  ) %>%
  mutate(log(Chao1),
         log(Shannon),
         log(Simpson),
         log(Observed)
         )
combined_indices_shrimp_trans <- combined_indices_shrimp_trans %>%
  mutate(`log(Shannon)` = `log(Shannon)` + 1.1703242,
         `log(Simpson)` = `log(Simpson)` + 2.3112746)

aof.shrimp <- function(x){
  y <- anova(aov(x ~ as.factor(day), data = combined_indices_shrimp_trans))
  return(y)
}

anova(aov(combined_indices_shrimp_trans$`log(Shannon)` ~ as.factor(combined_indices_shrimp_trans$treatment),
          data = combined_indices_shrimp_trans))
aof.plot <- function(x,y){
  plot <- ggplot(combined_indices_shrimp_trans, aes(x = treatment, y = x, fill = treatment)) +
    geom_boxplot() +
    # facet_wrap( ~ time, ncol = 2) +
    theme_classic() + 
    ylab(colnames(combined_indices_shrimp_trans[y])) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(plot)
}
# beta diversity
bc.shrimp <- as.matrix(distance(shrimp, method = "bray"))
pcoa_bc_shrimp <- ordinate(shrimp, "PCoA", "bray")
plot_ordination(shrimp,pcoa_bc_shrimp, color = "day", shape = "treatment")
# permutation ANOVA test
permanova_shrimp <- adonis2(bc.shrimp ~ day, data = shrimpexp, permutations = 999)
