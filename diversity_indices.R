# Load necessary libraries
library(phyloseq)
library(tidyverse)
library(SRS)

# Load the ASV and taxonomy data
asv_data <- read.csv("ASV.csv",row.names = 1)
asv_data <- asv_data[-106,] #remove negative control
taxonomy_data <- read.csv("taxa.csv",row.names = 1)
sample_data <- sample_data(read.csv("merge.csv",row.names = 1))
sample_data <- na.omit(sample_data)
time <- as.factor(sample_data$time.day.)
asv_data <- asv_data[match(rownames(sample_data),rownames(asv_data)),]
asv_data <- t(asv_data)
asv_data <- as.data.frame(asv_data)
asv_data_SRS <- SRS(asv_data,277) # what Cmin????
rownames(asv_data_SRS) <- rownames(asv_data)
asv_data_SRS <- t(asv_data_SRS)

# Convert taxonomy data to matrix and specify column names for Phyloseq
taxonomy_matrix <- as.matrix(taxonomy_data)
colnames(taxonomy_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# create phyloseq objects
otu_table <- otu_table(as.matrix(asv_data_SRS), taxa_are_rows = FALSE)
tax.table <- tax_table(taxonomy_matrix)
  
physeq <- phyloseq(otu_table,tax.table,sample_data)
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq) # ASV with not counts
t.genus <- tax_glom(physeq, taxrank = "Genus") # agglomerate genus
# find top 10 abundence genus
top10 <- names(sort(taxa_sums(t.genus),decreasing = TRUE))[1:10]
t.genus.prop <- transform_sample_counts(t.genus,function(x) x/sum(x))
t.genus.prop.top10 <- prune_taxa(top10,t.genus.prop)
## plot
t.genus.prop.top10@sam_data$time.day. <- factor(t.genus.prop.top10@sam_data$time.day.)
t.genus.prop.top10@sam_data <- t.genus.prop.top10@sam_data[order(t.genus.prop.top10@sam_data$time.day.)]
t.genus.prop.top10@otu_table <- t.genus.prop.top10@otu_table[match(rownames(t.genus.prop.top10@sam_data),rownames(t.genus.prop.top10@otu_table)),]

abudence_plot <- t.genus.prop.top10 %>%
  subset_samples(!rownames(t.genus.prop.top10@sam_data) %in% c("C2-0", "A3-3", "D1-6", "E3-6") &
    treatment %in% c("control")) %>% #"Nisin","Pediocin","Divergicin","MicrocinJ25"
  plot_bar(fill = "Genus") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_grid(time.day.~treatment, scales = "free") +
  geom_bar(aes(),stat = "identity",position = "stack") +
  theme_classic() +
  coord_flip()
  labs(title = "relative abudence of top 10 genus", 
       x = "samples", y = "relative abundence")

abudence_plot$data$Sample <- factor(abudence_plot$data$Sample,
                                    levels = rownames(sample_data),
                                    ordered = TRUE)
print(abudence_plot)
  dev.off()

# Calculate Chao1 diversity
bacterocin_exp <- subset_samples(physeq, 
                                 treatment %in% c("control","Nisin","Pediocin","Divergicin","MicrocinJ25"))

chao1_diversity <- estimate_richness(bacterocin_exp, measures = "Chao1")

# Add sample names as a column
chao1_diversity$Sample <- rownames(chao1_diversity)
chao1_diversity$treatment <- bacterocin_exp@sam_data$treatment

# Plot Chao1 diversity
chao1_plot <- ggplot(chao1_diversity, aes(x = Sample, y = Chao1, fill = treatment)) +
  geom_bar(stat = "identity")+
  theme_minimal() +
  xlab("Sample") +
  ylab("Chao1 Diversity Index") +
  coord_flip() +
  ggtitle("Chao1 Diversity Across Samples")

chao1_plot$data$Sample <- factor(chao1_plot$data$Sample,levels = chao1_plot$data$Sample,
                                 ordered = TRUE)
print(chao1_plot)
# Calculate Shannon diversity
shannon_diversity <- estimate_richness(bacterocin_exp, measures = "Shannon")

# Add sample names as a column
shannon_diversity$Sample <- rownames(shannon_diversity)
shannon_diversity$treatment <- bacterocin_exp@sam_data$treatment

# Plot Shannon diversity
shannon_plot <- ggplot(shannon_diversity, aes(x = Sample, y = Shannon, fill = treatment)) +
  coord_flip() +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Sample") +
  ylab("Shannon Diversity Index") +
  ggtitle("Shannon Diversity Across Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))# Adjust text angle for better readability if necessary

shannon_plot$data$Sample <- factor(shannon_plot$data$Sample,levels = shannon_plot$data$Sample,
                                   ordered = TRUE)
print(shannon_plot)
# Calculate Simpson diversity
simpson_diversity <- estimate_richness(bacterocin_exp, measures = "Simpson")

# Add sample names as a column
simpson_diversity$Sample <- rownames(simpson_diversity)
simpson_diversity$treatment <- bacterocin_exp@sam_data$treatment
# Plot Simpson diversity
simpson_diversity_plot <- 
  ggplot(simpson_diversity, aes(x = Sample, y = Simpson, fill = treatment)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  xlab("Sample") +
  ylab("Simpson Diversity Index") +
  ggtitle("Simpson Diversity Across Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

simpson_diversity_plot$data$Sample <-  factor(simpson_diversity_plot$data$Sample,
                                              levels = simpson_diversity_plot$data$Sample,
                                              ordered = TRUE)

print(simpson_diversity_plot)
# Adjust text angle for better readability if necessary

# Assuming 'physeq' is your phyloseq object
chao1_diversity <- estimate_richness(bacterocin_exp, measures = "Chao1")
shannon_diversity <- estimate_richness(bacterocin_exp, measures = "Shannon")
simpson_diversity <- estimate_richness(bacterocin_exp, measures = "Simpson")
observed_diversity <- estimate_richness(bacterocin_exp, measures = "Observed")

# Combine all indices into one data frame
combined_indices <- data.frame(Sample = rownames(chao1_diversity),
                               Chao1 = chao1_diversity$Chao1,
                               Shannon = shannon_diversity$Shannon,
                               Simpson = simpson_diversity$Simpson, 
                               Observed = observed_diversity$Observed)

# plots
plot_richness(bacterocin_exp, x = "treatment", measures = "Shannon") +
  geom_boxplot() +
  ylab("Shannon") +
  facet_grid(.~ time.day.)

plot_richness(bacterocin_exp, x = "treatment", measures = "Observed") +
  geom_boxplot() +
  ylab("Observed") +
  facet_grid(. ~ time.day.)

plot_richness(bacterocin_exp, x = "treatment", measures = "Simpson") +
  ylab("Simpson") +
  geom_boxplot() +
  facet_grid(. ~ time.day., scales = "free_y")

plot_richness(bacterocin_exp, x = "treatment", measures = "Chao1") +
  ylab("Chao1") +
  geom_boxplot() +
  facet_grid(. ~ time.day., scales = "free")
# View the head of the combined table
head(combined_indices)

# Export the combined table to a CSV file
write.csv(combined_indices, "combined_diversity_indices.csv", row.names = FALSE)

# beta diversity
bc <- distance(bacterocin_exp, method = "bray")

bc.tibble <- bc %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(time = as.numeric (str_extract(sample,"(?<=-)\\d+")),
         treatment = as.vector((str_extract(sample,"^[A-Za-z]")))) %>%
  # day_c = str_replace(sample, "c.*",""),
  # day_s = str_replace(sample, "s.*", "")) %>%
  group_by(time, sample, treatment) %>%
  summarize(median = median(value))%>%
  ungroup()

bc.tibble$replicate <- substr(bc.tibble$sample, 2,2)

bc.tibble$time <- sort(bc.tibble$time, decreasing = FALSE)

# add treatment label
bc.tibble <- bc.tibble %>%
  mutate(treatment = case_when(
    treatment == "A" ~ "control",
    treatment == "B" ~ "Nisin",
    treatment == "C" ~ "Pediocin",
    treatment == "D" ~ "Divergicin", 
    treatment == "E" ~ "MicrocinJ25"
  ))

## plot
ggplot(data = bc.tibble, aes(x=time, y=median, color = treatment,
                             group = paste0(replicate, treatment))) +
  geom_line(size = 0.25) +
  geom_smooth(aes(group = treatment), size = 3)


ordinate_bc <- ordinate(bacterocin_exp,method = "PCoA", distance = "bray")
plot_ordination(bacterocin_exp,ordinate_bc, color = "treatment") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = treatment))


sample_dat_bc <- data.frame(sample_data(bacterocin_exp))
sample_dat_bc$treatment <- factor(sample_dat_bc$treatment)
permeanova_result <- adonis2(bc ~ treatment + factor(time.day.), data = sample_dat_bc, permutations = 999)
