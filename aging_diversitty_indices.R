library(phyloseq)
library(SRS)
library(tidyverse)


asv_data <- read.csv("ASV.csv",row.names = 1)
asv_data <- asv_data[-106,] #remove negative control
taxonomy_data <- read.csv("taxa.csv",row.names = 1)
sample_data <- sample_data(read.csv("merge.csv",row.names = 1))
sample_data <- na.omit(sample_data)
time <- as.factor(sample_data$time.day.)
asv_data <- asv_data[match(rownames(sample_data),rownames(asv_data)),]
asv_data <- t(asv_data)
asv_data <- as.data.frame(asv_data)
asv_data_SRS <- SRS(asv_data,500) # what Cmin????
rownames(asv_data_SRS) <- rownames(asv_data)
asv_data_SRS <- t(asv_data_SRS)

taxonomy_matrix <- as.matrix(taxonomy_data)
colnames(taxonomy_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

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
  subset_samples(!rownames(t.genus.prop.top10@sam_data) %in% c("oc2a", "12c1a", "18c2c", "6c3c", "9c2a") &
                   treatment %in% c("c","s")) %>% #"Nisin","Pediocin","Divergicin","MicrocinJ25"
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


#### alpha diversity
aging_experiment <- subset_samples(physeq, 
                                   treatment %in% c("c","s"))

chao1_diversity <- estimate_richness(aging_experiment, measures = "Chao1")
shannon_diversity <- estimate_richness(aging_experiment, measures = "Shannon")
simpson_diversity <- estimate_richness(aging_experiment, measures = "Simpson")
observed_diversity <- estimate_richness(aging_experiment, measures = "Observed")

combined_indices <- data.frame(Sample = rownames(chao1_diversity),
                               Chao1 = chao1_diversity$Chao1,
                               Shannon = shannon_diversity$Shannon,
                               Simpson = simpson_diversity$Simpson, 
                               Observed = observed_diversity$Observed)

plot_richness(aging_experiment, x = "treatment", measures = "Shannon") +
  geom_boxplot() +
  ylab("Shannon") +
  facet_grid(.~ time.day.)

plot_richness(aging_experiment, x = "treatment", measures = "Observed") +
  geom_boxplot() +
  ylab("Observed") +
  facet_grid(. ~ time.day.)

plot_richness(aging_experiment, x = "treatment", measures = "Simpson") +
  ylab("Simpson") +
  geom_boxplot() +
  facet_grid(. ~ time.day., scales = "free_y")

plot_richness(aging_experiment, x = "treatment", measures = "Chao1") +
  ylab("Chao1") +
  geom_boxplot() +
  facet_grid(. ~ time.day., scales = "free")

write.csv(combined_indices, "aging_combined_diversity_indices.csv", row.names = FALSE)


bc <- distance(aging_experiment, method = "bray")

bc.tibble <- bc %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(coho = str_replace(sample, ".*c.*", "coho"),
         sockeye = str_replace(sample, ".*s.*","sockeye"),
         time = as.numeric (str_replace(sample,"c.*|s.*", "")),
         treatment = as.vector(str_match(sample,"[cs]"))) %>%
         # day_c = str_replace(sample, "c.*",""),
         # day_s = str_replace(sample, "s.*", "")) %>%
  group_by(time, sample, treatment) %>%
  summarize(median = median(value))%>%
  ungroup()

bc.tibble$replicate <- substr(bc.tibble$sample, nchar(bc.tibble$sample) -1, 
                              nchar(bc.tibble$sample))
bc.tibble$replicate[19] <- "1a"

bc.tibble$time <- sort(bc.tibble$time, decreasing = FALSE)

ggplot(data = bc.tibble, aes(x=time, y=median, color = treatment,
                             group = paste0(replicate, treatment))) +
  geom_line(size = 0.25) +
  geom_smooth(aes(group = treatment), size = 3)


ordinate_bc <- ordinate(aging_experiment,method = "PCoA", distance = "bray")
plot_ordination(aging_experiment,ordinate_bc, color = "time.day.", shape = "treatment") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = treatment))

plot_ordination(aging_experiment,ordinate_bc, color = "treatment") +
  geom_point(size = 4) +
  stat_ellipse(aes(group = treatment))


sample_dat_bc <- data.frame(sample_data(aging_experiment))
sample_dat_bc$treatment <- factor(sample_dat_bc$treatment)
permeanova_result <- adonis2(bc ~ treatment + factor(time.day.), data = sample_dat_bc, permutations = 999)
