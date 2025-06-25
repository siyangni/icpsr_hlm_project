#####
## Principal Component Analysis for Delinquency Variables
## ICP-H Project
#####

library(pacman)
p_load(dplyr, readr, stats, psych, factoextra, corrplot, car, ggplot2, gridExtra)

#####
## Load and prepare data
#####

# Load processed dataset
df <- read_csv("data/processed.csv")

# Recode B24 variables from 1,2,3,4 to 0,1,2,3  
df <- df %>%
  mutate(
    B24j_recoded = case_when(
      B24j == 1 ~ 0,
      B24j == 2 ~ 1,
      B24j == 3 ~ 2,
      B24j == 4 ~ 3,
      TRUE ~ NA_real_
    ),
    B24k_recoded = case_when(
      B24k == 1 ~ 0,
      B24k == 2 ~ 1,
      B24k == 3 ~ 2,
      B24k == 4 ~ 3,
      TRUE ~ NA_real_
    ),
    B24l_recoded = case_when(
      B24l == 1 ~ 0,
      B24l == 2 ~ 1,
      B24l == 3 ~ 2,
      B24l == 4 ~ 3,
      TRUE ~ NA_real_
    ),
    B24m_recoded = case_when(
      B24m == 1 ~ 0,
      B24m == 2 ~ 1,
      B24m == 3 ~ 2,
      B24m == 4 ~ 3,
      TRUE ~ NA_real_
    ),
    B24n_recoded = case_when(
      B24n == 1 ~ 0,
      B24n == 2 ~ 1,
      B24n == 3 ~ 2,
      B24n == 4 ~ 3,
      TRUE ~ NA_real_
    )
  )

# Define delinquency variables (from eda.R lines 63-65)
delinq.vars <- c("B24j_recoded", "B24k_recoded", "B24l_recoded", "B24m_recoded", "B24n_recoded", 
                 "B27a", "B27b", "B27c", "B27d", "B27e", 
                 "B27f", "B27g", "B27h", "B27i", "B27j")

# Prepare data for PCA
df_pca <- df[, delinq.vars]
df_pca <- na.omit(df_pca)

cat("=== DATA PREPARATION ===\n")
cat("Sample size for PCA:", nrow(df_pca), "\n")
cat("Number of variables:", ncol(df_pca), "\n\n")

# Check data structure before standardization
cat("=== VARIABLE SUMMARIES (Before Standardization) ===\n")
summary(df_pca)

# Standardize all variables (mean = 0, sd = 1)
cat("\n=== STANDARDIZING VARIABLES ===\n")
df_pca_raw <- df_pca  # Keep original for reference
df_pca_scaled <- scale(df_pca, center = TRUE, scale = TRUE)
df_pca <- as.data.frame(df_pca_scaled)

cat("Variables standardized (centered and scaled)\n")
cat("All variables now have mean ≈ 0 and standard deviation = 1\n\n")

# Check data structure after standardization
cat("=== VARIABLE SUMMARIES (After Standardization) ===\n")
summary(df_pca)

# Verify standardization
cat("\nStandardization verification:\n")
cat("Means (should be close to 0):\n")
print(round(colMeans(df_pca), 6))
cat("\nStandard deviations (should be 1):\n")
print(round(apply(df_pca, 2, sd), 6))

#####
## Preliminary checks for PCA suitability
#####

cat("\n=== PCA SUITABILITY TESTS ===\n")

# Correlation matrix
cor_matrix <- cor(df_pca)
cat("Correlation matrix computed\n")

# Bartlett's test of sphericity
bartlett_result <- cortest.bartlett(cor_matrix, n = nrow(df_pca))
cat("Bartlett's Test of Sphericity:\n")
cat("Chi-square =", bartlett_result$chisq, "\n")
cat("p-value =", bartlett_result$p.value, "\n")
cat("Interpretation:", ifelse(bartlett_result$p.value < 0.05, 
                             "Suitable for PCA (p < 0.05)", 
                             "Not suitable for PCA (p >= 0.05)"), "\n\n")

# Kaiser-Meyer-Olkin (KMO) test
kmo_result <- KMO(df_pca)
cat("Kaiser-Meyer-Olkin (KMO) Test:\n")
cat("Overall MSA =", kmo_result$MSA, "\n")
cat("Interpretation:", 
    ifelse(kmo_result$MSA >= 0.8, "Excellent", 
           ifelse(kmo_result$MSA >= 0.7, "Good",
                  ifelse(kmo_result$MSA >= 0.6, "Mediocre",
                         ifelse(kmo_result$MSA >= 0.5, "Poor", "Unacceptable")))), "\n")

cat("Individual MSA values:\n")
print(round(kmo_result$MSAi, 3))

#####
## Principal Component Analysis
#####

cat("\n=== PRINCIPAL COMPONENT ANALYSIS ===\n")

# Perform PCA using prcomp (variables already standardized)
pca_result <- prcomp(df_pca, center = FALSE, scale. = FALSE)

# Summary of PCA
cat("PCA Summary:\n")
summary(pca_result)

# Extract eigenvalues
eigenvalues <- pca_result$sdev^2
cat("\nEigenvalues:\n")
print(round(eigenvalues, 3))

# Proportion of variance explained
prop_var <- eigenvalues / sum(eigenvalues)
cumulative_var <- cumsum(prop_var)

cat("\nProportion of variance explained:\n")
variance_table <- data.frame(
  Component = 1:length(eigenvalues),
  Eigenvalue = round(eigenvalues, 3),
  Prop_Variance = round(prop_var, 3),
  Cumulative_Variance = round(cumulative_var, 3)
)
print(variance_table)

#####
## Component retention criteria
#####

cat("\n=== COMPONENT RETENTION CRITERIA ===\n")

# Kaiser criterion (eigenvalue > 1)
kaiser_components <- sum(eigenvalues > 1)
cat("Kaiser criterion (eigenvalue > 1):", kaiser_components, "components\n")

# Scree test (visual inspection needed)
cat("Scree test: See scree plot for visual inspection\n")

# Variance explained criteria
var_80 <- sum(cumulative_var <= 0.8) + 1
var_90 <- sum(cumulative_var <= 0.9) + 1
cat("Components explaining 80% variance:", var_80, "\n")
cat("Components explaining 90% variance:", var_90, "\n")

# Parallel analysis
parallel_result <- fa.parallel(df_pca, fa = "pc", n.iter = 100, show.legend = FALSE)
parallel_components <- parallel_result$ncomp
cat("Parallel analysis suggests:", parallel_components, "components\n")

#####
## Visualizations
#####

cat("\n=== CREATING VISUALIZATIONS ===\n")

# 1. Scree plot
scree_plot <- fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50),
                       title = "Scree Plot - Delinquency Variables PCA")
print(scree_plot)

# 2. Correlation matrix heatmap
corrplot(cor_matrix, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.8, tl.col = "black")
title("Correlation Matrix - Delinquency Variables")

# 3. Biplot - first two components
biplot_12 <- fviz_pca_biplot(pca_result, 
                             title = "PCA Biplot - PC1 vs PC2",
                             geom.ind = "point",
                             pointsize = 0.5,
                             alpha.ind = 0.3,
                             col.var = "contrib",
                             gradient.cols = c("white", "blue", "red"))
print(biplot_12)

# 4. Variable contributions to PC1 and PC2
contrib_pc1 <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 15,
                            title = "Contribution to PC1")
contrib_pc2 <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 15,
                            title = "Contribution to PC2")
print(contrib_pc1)
print(contrib_pc2)

# 5. Quality of representation (cos2)
cos2_plot <- fviz_pca_var(pca_result, col.var = "cos2",
                          gradient.cols = c("black", "orange", "red"),
                          repel = TRUE,
                          title = "Variables - Quality of Representation")
print(cos2_plot)

#####
## Component loadings analysis
#####

cat("\n=== COMPONENT LOADINGS ANALYSIS ===\n")

# Extract loadings for first few components
n_components <- min(5, ncol(pca_result$rotation))
loadings_matrix <- pca_result$rotation[, 1:n_components]

cat("Component loadings (first", n_components, "components):\n")
print(round(loadings_matrix, 3))

# Identify variables with high loadings (>0.3 or <-0.3) for each component
cat("\nHigh loadings (|loading| > 0.3) by component:\n")
for(i in 1:n_components) {
  high_loadings <- which(abs(loadings_matrix[, i]) > 0.3)
  if(length(high_loadings) > 0) {
    cat(paste("PC", i, ":\n", sep = ""))
    for(j in high_loadings) {
      cat("  ", rownames(loadings_matrix)[j], ":", 
          round(loadings_matrix[j, i], 3), "\n")
    }
  }
}

#####
## Component scores
#####

cat("\n=== COMPONENT SCORES ===\n")

# Extract component scores
component_scores <- pca_result$x[, 1:n_components]

# Add component scores to original dataframe indices
df_with_scores <- df_pca
for(i in 1:n_components) {
  df_with_scores[paste0("PC", i)] <- component_scores[, i]
}

cat("Component scores calculated and added to dataset\n")
cat("First few component scores:\n")
print(head(component_scores))

# Summary statistics for component scores
cat("\nSummary statistics for component scores:\n")
print(summary(component_scores))

#####
## Interpretation and recommendations
#####

cat("\n=== INTERPRETATION AND RECOMMENDATIONS ===\n")

cat("Sample size:", nrow(df_pca), "\n")
cat("Number of variables:", length(delinq.vars), "\n")
cat("KMO Overall MSA:", round(kmo_result$MSA, 3), "\n")
cat("Bartlett's test p-value:", format(bartlett_result$p.value, scientific = TRUE), "\n")
cat("Components with eigenvalue > 1:", kaiser_components, "\n")
cat("Components suggested by parallel analysis:", parallel_components, "\n")
cat("Variance explained by first component:", round(prop_var[1] * 100, 1), "%\n")
cat("Cumulative variance explained by first", kaiser_components, "components:", 
    round(cumulative_var[kaiser_components] * 100, 1), "%\n")

# Recommendations
cat("\nRecommendations:\n")
if(kmo_result$MSA >= 0.6 & bartlett_result$p.value < 0.05) {
  cat("- Data is suitable for PCA based on KMO and Bartlett's tests\n")
} else {
  cat("- Data may not be well-suited for PCA - consider alternative approaches\n")
}

if(kaiser_components == parallel_components) {
  cat("- Kaiser criterion and parallel analysis agree on", kaiser_components, "components\n")
} else {
  cat("- Kaiser criterion suggests", kaiser_components, "components\n")
  cat("- Parallel analysis suggests", parallel_components, "components\n")
  cat("- Consider parallel analysis result as more reliable\n")
}

cat("- First component explains", round(prop_var[1] * 100, 1), 
    "% of variance - indicates", 
    ifelse(prop_var[1] > 0.4, "strong general factor", "moderate general factor"), "\n")

#####
## Save results
#####

cat("\n=== SAVING RESULTS ===\n")

# Save PCA object
saveRDS(pca_result, "pca_delinquency_result.rds")

# Save component scores
write_csv(df_with_scores, "data/processed_with_pca_scores.csv")

# Save loadings matrix
write.csv(loadings_matrix, "pca_loadings_matrix.csv")

# Save variance explained table
write.csv(variance_table, "pca_variance_explained.csv")

cat("Results saved:\n")
cat("- PCA object: pca_delinquency_result.rds\n")
cat("- Data with PCA scores: data/processed_with_pca_scores.csv\n")
cat("- Loadings matrix: pca_loadings_matrix.csv\n")
cat("- Variance table: pca_variance_explained.csv\n")

cat("\nPCA analysis complete!\n") 