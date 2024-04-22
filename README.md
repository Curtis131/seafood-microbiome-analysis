# BINF531 term project

**before running the code, please:**
- make sure you've installed biocondutor;
 ``if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")``

- make sure you've installed DADA2:
`` BiocManager::install("dada2")``

- make sure you've changed the working directory on the 1st line.

- make sure you've placed the [results](https://drive.google.com/drive/folders/1ZSMh-QzIVVQmPsmO3HnX2cv8G53r1yc8?usp=drive_link) in the working directory

## Description for files

> ``term_project.R`` contains preprocessing pipeline using DADA2.

> ``preprocessingresult.csv`` contains preprocessing result.

> ``create_metadata.R`` contains script used to generate experiment metadata.

> ``diversity_indices.R`` contains script used for analysis & visualization of alpha & beta diversity for the bacterocin experiment.

> ``diversity_indices_analysis.R`` contains script used for probability testing of diversity indices.
