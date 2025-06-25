#####
## Exploratory Data Analysis for ICP-H Project
## Delinquency Scale Development using Mokken Scale Analysis
#####

library(pacman)
p_load(stats, car, lattice, foreign, mokken, mirt, dplyr, readr, ltm, magrittr)

#####
## Load processed dataset
#####

df <- read_csv("data/processed.csv")
head(df)

#####
## Recode delinquency variables for proper ordinal analysis
## B24j, B24k, B24l, B24m, B24n: recode 1,2,3,4 to 3,2,1,0
## B27a-B27j: use original variables as is
#####

# Recode B24 variables to reverse the scale (1,2,3,4 -> 0, 1, 2, 3) 
df <- df %>%
  mutate(
    B24j_recoded = case_when(
      B24j == 1 ~ 3,
      B24j == 2 ~ 2,
      B24j == 3 ~ 1,
      B24j == 4 ~ 0,
      TRUE ~ NA_real_
    ),
    B24k_recoded = case_when(
      B24k == 1 ~ 3,
      B24k == 2 ~ 2,
      B24k == 3 ~ 1,
      B24k == 4 ~ 0,
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

# Define delinquency variables for analysis
delinq.vars <- c("B24j_recoded", "B24k_recoded", "B24l_recoded", "B24m_recoded", "B24n_recoded", 
                 "B27a", "B27b", "B27c", "B27d", "B27e", 
                 "B27f", "B27g", "B27h", "B27i", "B27j")

# Summary of all delinquency variables
cat("=== DELINQUENCY VARIABLES SUMMARY ===\n")
cat("B24 variables have been recoded: 1,2,3,4 -> 3,2,1,0\n")
cat("B27 variables are used in original form\n\n")
for(var in delinq.vars) {
  if(var %in% names(df)) {
    cat("Variable:", var, "\n")
    print(table(df[[var]], useNA = "ifany"))
    cat("\n")
  }
}

#####
## Subset data and remove missing values
## for delinquency scale analysis
#####

# Select variables using base R approach for compatibility
df_clean <- df[, delinq.vars]
df_clean <- na.omit(df_clean)

cat("Sample size after removing missing data:", nrow(df_clean), "\n")

#####
## Mokken Scale Analysis for Delinquency Variables
#####

# Step 1: Automated Item Selection Procedure (AISP)
cat("=== AUTOMATED ITEM SELECTION PROCEDURE ===\n")
aisp(df_clean[,delinq.vars], verbose = TRUE)

# Step 2: Check scalability coefficients
cat("\n=== SCALABILITY COEFFICIENTS ===\n")
coefH(df_clean[,delinq.vars])

# Step 3: Check assumptions of Mokken scaling

# Check monotonicity assumption
cat("\n=== MONOTONICITY CHECK ===\n")
mono_check <- check.monotonicity(df_clean[,delinq.vars])
summary(mono_check)

# Check restscore assumption  
cat("\n=== RESTSCORE CHECK ===\n")
rest_check <- check.restscore(df_clean[,delinq.vars])
summary(rest_check)

# Check P-matrix (pairwise scalability)
cat("\n=== P-MATRIX CHECK ===\n")
pmat_check <- check.pmatrix(df_clean[,delinq.vars])
summary(pmat_check)

# Check Invariant Item Ordering (IIO)
cat("\n=== INVARIANT ITEM ORDERING CHECK ===\n")
iio_check <- check.iio(df_clean[,delinq.vars])
summary(iio_check)

# Check reliability
cat("\n=== RELIABILITY CHECK ===\n")
reliability <- check.reliability(df_clean[,delinq.vars])
print(reliability)

# Create delinquency scale score (sum of items)
df$delinquency_mokken <- rowSums(df[,delinq.vars], na.rm = TRUE)

#####
## If any items are problematic, create revised scale
## Check if we need to remove any items based on AISP results
#####

# Assuming we might need to remove some items based on violations
# (This will depend on your actual results - you may need to adjust)

# Example: If certain items violate assumptions, create a revised scale
# delinq.vars.revised <- c("b24j", "b24k", "b24l", "b24m", "b24n", 
#                          "b27a", "b27b", "b27c", "b27d", "b27e")

# Uncomment and modify based on your results:
# cat("\n=== REVISED SCALE ANALYSIS ===\n")
# coefH(df_clean[,delinq.vars.revised])
# summary(check.monotonicity(df_clean[,delinq.vars.revised]))
# summary(check.restscore(df_clean[,delinq.vars.revised]))
# summary(check.pmatrix(df_clean[,delinq.vars.revised]))
# summary(check.iio(df_clean[,delinq.vars.revised]))
# check.reliability(df_clean[,delinq.vars.revised])

#####
## Visualizations
#####

# Distribution of delinquency scale
histogram(~as.factor(delinquency_mokken),
          data = df,
          xlab = "Delinquency Scale Score",
          main = "Distribution of Delinquency Scale"
)

# Item difficulty ranking (proportion endorsing each item)
item_props <- colMeans(df_clean[,delinq.vars], na.rm = TRUE)
cat("\n=== ITEM DIFFICULTY RANKING ===\n")
print(sort(item_props, decreasing = TRUE))

#####
## Parametric IRT Analysis for comparison
#####

cat("\n=== PARAMETRIC IRT ANALYSIS ===\n")

# Empirical plots
empirical_plot(df_clean[,delinq.vars], smooth = TRUE,
               which.items = 1:length(delinq.vars))

# Rasch model
cat("\n--- Rasch Model ---\n")
rasch.delinq <- mirt(df_clean[,delinq.vars], 1, itemtype = "Rasch", SE = TRUE)
print(M2(rasch.delinq))
itemfit_rasch <- itemfit(rasch.delinq)
print(itemfit_rasch)

# Plot item characteristic curves for Rasch model
plot(rasch.delinq, type = "trace", which.items = 1:length(delinq.vars),
     main = "Item Characteristic Curves - Rasch Model")

# 2PL model
cat("\n--- 2PL Model ---\n")
twopl.delinq <- mirt(df_clean[,delinq.vars], 1, itemtype = "2PL", SE = TRUE)
print(M2(twopl.delinq))
itemfit_2pl <- itemfit(twopl.delinq)
print(itemfit_2pl)

# Plot item characteristic curves for 2PL model
plot(twopl.delinq, type = "trace", which.items = 1:length(delinq.vars),
     main = "Item Characteristic Curves - 2PL Model")

# Get IRT parameters
cat("\n--- 2PL Parameters ---\n")
irt_params <- coef(twopl.delinq, IRTpars = TRUE)
print(irt_params)

# 3PL model
cat("\n--- 3PL Model ---\n")
tryCatch({
  threepl.delinq <- mirt(df_clean[,delinq.vars], 1, itemtype = "3PL")
  print(M2(threepl.delinq))
  print(coef(threepl.delinq, IRTpars = TRUE))
  print(itemfit(threepl.delinq))
}, error = function(e) {
  cat("3PL model failed to converge or had issues:", e$message, "\n")
})

# Model comparison
cat("\n--- Model Comparison ---\n")
print(anova(rasch.delinq, twopl.delinq))

#####
## Retrieve and compare scores from different methods
#####

# Get theta scores from 2PL model
theta_scores_full <- fscores(twopl.delinq, full.scores.SE = TRUE)[,1]

# Create a temporary dataframe to handle dimension mismatch
df_with_theta <- df_clean %>%
  mutate(theta_scores = theta_scores_full,
         delinquency_mokken = rowSums(df_clean[,delinq.vars], na.rm = TRUE))

# Compare Mokken sum scores with IRT theta scores
correlation <- cor(df_with_theta$delinquency_mokken, df_with_theta$theta_scores, use = "complete.obs")
cat("\nCorrelation between Mokken sum scores and 2PL theta scores:", correlation, "\n")

# Scatterplot comparing scores
set.seed(1234)
xyplot(jitter(theta_scores) ~ jitter(delinquency_mokken),
       data = df_with_theta,
       aspect = 1,
       col = "black",
       xlab = "Mokken Sum Scores",
       ylab = "2PL Theta Scores",
       main = "Comparison of Mokken and IRT Scores"
)

#####
## Calculate point-biserial correlations
## (similar to discrimination parameters)
#####

cat("\n=== POINT-BISERIAL CORRELATIONS ===\n")

# Calculate point-biserial correlations for each item with total scale
for(i in 1:length(delinq.vars)) {
  item_name <- delinq.vars[i]
  # Calculate scale without current item
  scale_without_item <- rowSums(df_clean[,-i], na.rm = TRUE)
  
  # Calculate point-biserial correlation
  pb_cor <- biserial.cor(scale_without_item, df_clean[,i], level = 2)
  cat(item_name, ":", round(pb_cor, 3), "\n")
}

#####
## Summary and recommendations
#####

cat("\n=== SUMMARY ===\n")
cat("Number of items in delinquency scale:", length(delinq.vars), "\n")
cat("Sample size for analysis:", nrow(df_clean), "\n")
cat("Correlation between Mokken and IRT scores:", round(correlation, 3), "\n")

# Save results for further analysis
write_csv(df, "data/processed_with_scales.csv")

cat("\nAnalysis complete. Results saved to data/processed_with_scales.csv\n")
