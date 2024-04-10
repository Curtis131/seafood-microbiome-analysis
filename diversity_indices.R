# Load necessary libraries
library(phyloseq)
library(ggplot2)

# Load the ASV and taxonomy data
asv_data <- read.csv("seqtab.nochim.csv", row.names = 1)
taxonomy_data <- read.csv("taxa.csv", row.names = 1)

# Convert taxonomy data to matrix and specify column names for Phyloseq
taxonomy_matrix <- as.matrix(taxonomy_data)
colnames(taxonomy_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Create an otu_table object from the ASV data
otu_table <- otu_table(as.matrix(asv_data), taxa_are_rows = FALSE)

# Create a tax_table object from the taxonomy data
tax_table <- tax_table(taxonomy_matrix)

# Create the Phyloseq object
physeq <- phyloseq(otu_table, tax_table)

# Calculate Chao1 diversity
chao1_diversity <- estimate_richness(physeq, measures = "Chao1")

# Add sample names as a column
chao1_diversity$Sample <- rownames(chao1_diversity)

# Plot Chao1 diversity
ggplot(chao1_diversity, aes(x = Sample, y = Chao1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Sample") +
  ylab("Chao1 Diversity Index") +
  ggtitle("Chao1 Diversity Across Samples")

# Calculate Shannon diversity
shannon_diversity <- estimate_richness(physeq, measures = "Shannon")

# Add sample names as a column
shannon_diversity$Sample <- rownames(shannon_diversity)

# Plot Shannon diversity
ggplot(shannon_diversity, aes(x = Sample, y = Shannon)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Sample") +
  ylab("Shannon Diversity Index") +
  ggtitle("Shannon Diversity Across Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for better readability if necessary

# Calculate Simpson diversity
simpson_diversity <- estimate_richness(physeq, measures = "Simpson")

# Add sample names as a column
simpson_diversity$Sample <- rownames(simpson_diversity)
# Plot Simpson diversity
ggplot(simpson_diversity, aes(x = Sample, y = Simpson)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Sample") +
  ylab("Simpson Diversity Index") +
  ggtitle("Simpson Diversity Across Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angle for better readability if necessary

# Assuming 'physeq' is your phyloseq object
chao1_diversity <- estimate_richness(physeq, measures = "Chao1")
shannon_diversity <- estimate_richness(physeq, measures = "Shannon")
simpson_diversity <- estimate_richness(physeq, measures = "Simpson")

# Combine all indices into one data frame
combined_indices <- data.frame(Sample = rownames(chao1_diversity),
                               Chao1 = chao1_diversity$Chao1,
                               Shannon = shannon_diversity$Shannon,
                               Simpson = simpson_diversity$Simpson)
# View the head of the combined table
head(combined_indices)

# Export the combined table to a CSV file
write.csv(combined_indices, "combined_diversity_indices.csv", row.names = FALSE)

