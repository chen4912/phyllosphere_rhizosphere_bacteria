###################Random Sampling###################
set.seed(123123)
library(vegan)
library(dplyr)
ASV <- read.xlsx("field_survey.xlsx", sheet = "ASV", colNames = TRUE, rowNames = TRUE)
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

species_list <- c("F", "P", "T")
compartments <- c("L", "R")
months <- c("Sep", "Oct", "Dec", "Jun", "Aug")

results <- data.frame(Species = character(),
                      Compartment = character(),
                      Month = character(),
                      Significant_Count = integer(),
                      Non_Significant_Count = integer(),
                      stringsAsFactors = FALSE)

for (Species in species_list) {
  for (Compartment in compartments) {
    for (Month in months) {
      compartment_data <- subset(group, species == Species & compartment == Compartment & month == Month)
      filtered_data <- ASV[, colnames(ASV) %in% rownames(compartment_data)]
      otu_data <- filtered_data[rowSums(filtered_data > 0) > 0, ] 
      
      original_community <- as.matrix(otu_data)
      dist_original <- as.matrix(vegdist(t(original_community), method = "bray")) 
      
      n_samples <- 1000
      significant_results <- numeric(n_samples)
      
      n_otu <- ncol(otu_data)
      for (i in 1:n_samples) {
        sampled_otu_indices <- sample(1:n_otu, size = n_otu, replace = TRUE)
        sampled_otu <- otu_data[, sampled_otu_indices]
        random_otu <- merge(otu_data, sampled_otu, by = "row.names", all = TRUE)
        rownames(random_otu) <- random_otu[,1]
        random_otu <- random_otu[,-1]
        random_otu[is.na(random_otu)] <- 0
        dist_random <- vegdist(t(random_otu), method = "bray")
        group1 <- factor(c(rep("Random", ncol(sampled_otu)), rep("Original", ncol(otu_data))))
        # Perform PERMANOVA analysis
        permanova_result <- adonis2(dist_random ~ group1)
        permanova_result$`Pr(>F)`
        significant_results[i] <- ifelse(permanova_result$`Pr(>F)`[1] < 0.05, 1, 0)
      }
      
      significant_count <- sum(significant_results)
      non_significant_count <- n_samples - significant_count
      
      results <- rbind(results, data.frame(Species = Species,
                                           Compartment = Compartment,
                                           Month = Month,
                                           Significant_Count = significant_count,
                                           Non_Significant_Count = non_significant_count,
                                           stringsAsFactors = FALSE))
    }
  }
}

print(results)
