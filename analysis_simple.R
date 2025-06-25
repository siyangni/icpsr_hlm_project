# ===============================================================================
# Three-Level Hierarchical Negative Binomial Regression Analysis
# ICP-H Project - Simplified Version
# ===============================================================================

# Load only essential packages
library(dplyr)
library(readr)
library(glmmTMB)
library(broom.mixed)

# ===============================================================================
# LOAD AND PREPARE DATA
# ===============================================================================

# Load processed data
df <- read_csv("data/processed.csv")

# Check data structure
cat("=== DATA OVERVIEW ===\n")
cat("Sample size:", nrow(df), "\n")
cat("Number of areas:", length(unique(df$area)), "\n")
cat("Number of schools:", length(unique(df$school)), "\n\n")

# Define model variables
model_vars <- c("delinquency", "mgb", "stegb", "ltegb", "vbs", "nrwa", "nrws", 
                "nrwp", "pf", "np", "nle", "sne", "apc", "ats", "atp", "D45", 
                "D47", "bcv", "ica", "delinf", "fses", "age", "gender", "lbc1",
                "school", "area")

# Create analysis dataset with complete cases
df_model <- df %>%
  select(all_of(model_vars)) %>%
  na.omit() %>%
  mutate(
    gender = factor(gender, levels = c(0, 1), labels = c("Male", "Female")),
    lbc1 = factor(lbc1, levels = c(0, 1, 2, 3), 
                  labels = c("Non-LBC", "Father-migrant", "Mother-migrant", "Both-migrant")),
    school = factor(school),
    area = factor(area)
  )

cat("Analysis sample size:", nrow(df_model), "\n")
cat("Percentage retained:", round(100 * nrow(df_model) / nrow(df), 2), "%\n\n")

# ===============================================================================
# MODEL SPECIFICATION
# ===============================================================================

# Base formula from ch4_regression.R lines 103-107
formula_base <- "delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
                 nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
                 fses + age + gender + lbc1"

# Add interaction terms with lbc1
interactions <- c("lbc1:mgb", "lbc1:stegb", "lbc1:ltegb", "lbc1:vbs", "lbc1:nrwa", 
                  "lbc1:nrwp", "lbc1:pf", "lbc1:np", "lbc1:nle")

formula_with_interactions <- paste(formula_base, "+", paste(interactions, collapse = " + "))

# Add three-level random effects structure
formula_full <- paste(formula_with_interactions, "+ (1|area/school)")

cat("=== MODEL FORMULA ===\n")
cat(formula_full, "\n\n")

# ===============================================================================
# FIT MODELS
# ===============================================================================

cat("=== FITTING MODELS ===\n")

# Single-level model for comparison
cat("Fitting single-level model...\n")
model_1level <- glmmTMB(
  formula = as.formula(formula_with_interactions),
  data = df_model,
  family = nbinom2()
)

# Three-level hierarchical model
cat("Fitting three-level model...\n")
model_3level <- glmmTMB(
  formula = as.formula(formula_full),
  data = df_model,
  family = nbinom2()
)

# ===============================================================================
# MODEL RESULTS
# ===============================================================================

cat("\n=== MODEL COMPARISON ===\n")
cat("Single-level AIC:", AIC(model_1level), "\n")
cat("Three-level AIC:", AIC(model_3level), "\n")
cat("AIC improvement:", AIC(model_1level) - AIC(model_3level), "\n\n")

# Model summary
cat("=== THREE-LEVEL MODEL SUMMARY ===\n")
summary(model_3level)

# ===============================================================================
# EXTRACT RESULTS
# ===============================================================================

# Fixed effects
fixed_effects <- broom.mixed::tidy(model_3level, effects = "fixed")
cat("\n=== FIXED EFFECTS ===\n")
print(fixed_effects)

# Calculate IRRs
fixed_effects$IRR <- exp(fixed_effects$estimate)
fixed_effects$IRR_lower <- exp(fixed_effects$estimate - 1.96 * fixed_effects$std.error)
fixed_effects$IRR_upper <- exp(fixed_effects$estimate + 1.96 * fixed_effects$std.error)

cat("\n=== INCIDENT RATE RATIOS ===\n")
print(fixed_effects %>% 
      select(term, IRR, IRR_lower, IRR_upper, p.value) %>%
      mutate(across(c(IRR, IRR_lower, IRR_upper), ~round(.x, 3))))

# Random effects
random_effects <- broom.mixed::tidy(model_3level, effects = "ran_pars")
cat("\n=== RANDOM EFFECTS ===\n")
print(random_effects)

# Calculate ICCs
var_area <- random_effects$estimate[random_effects$term == "sd__(Intercept).area"]^2
var_school <- random_effects$estimate[random_effects$term == "sd__(Intercept).school:area"]^2
var_residual <- pi^2/3  # Approximation for negative binomial

total_var <- var_area + var_school + var_residual
icc_area <- var_area / total_var
icc_school <- var_school / total_var

cat("\n=== INTRACLASS CORRELATIONS ===\n")
cat("ICC Area (Level 3):", round(icc_area, 4), "\n")
cat("ICC School (Level 2):", round(icc_school, 4), "\n")
cat("ICC Combined:", round((var_area + var_school) / total_var, 4), "\n")

# ===============================================================================
# INTERACTION EFFECTS
# ===============================================================================

# Extract interaction terms
interaction_effects <- fixed_effects %>%
  filter(grepl("lbc1.*:", term)) %>%
  arrange(p.value)

cat("\n=== INTERACTION EFFECTS ===\n")
cat("Significant interactions (p < 0.05):\n")
significant_interactions <- interaction_effects %>% filter(p.value < 0.05)

if(nrow(significant_interactions) > 0) {
  print(significant_interactions %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value))
} else {
  cat("No significant interactions found.\n")
}

cat("\nAll interaction effects:\n")
print(interaction_effects %>% 
      select(term, IRR, IRR_lower, IRR_upper, p.value))

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

# Save model objects
saveRDS(model_3level, "model_3level_nb.rds")
saveRDS(model_1level, "model_1level_nb.rds")

# Save results tables
write_csv(fixed_effects, "fixed_effects_results.csv")
write_csv(random_effects, "random_effects_results.csv")

cat("\n=== RESULTS SAVED ===\n")
cat("- Three-level model: model_3level_nb.rds\n")
cat("- Single-level model: model_1level_nb.rds\n")
cat("- Fixed effects: fixed_effects_results.csv\n")
cat("- Random effects: random_effects_results.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n") 