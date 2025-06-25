# ===============================================================================
# Three-Level Hierarchical Negative Binomial Regression Analysis (Simplified)
# ICP-H Project - Using lbc0 instead of lbc1
# 
# Level 1: Individual (main effects only, no interactions)
# Level 2: School (random intercept only)
# Level 3: Area (random intercept only)
# ===============================================================================

# Load required packages
library(pacman)
p_load(dplyr, readr, lme4, glmmTMB, sjPlot, sjstats, performance, 
       car, broom.mixed, dotwhisker)

# ===============================================================================
# LOAD DATA
# ===============================================================================

# Load processed data from data_preprocessing.R
df <- read_csv("data/processed.csv")

# Check data structure
cat("=== DATA OVERVIEW ===\n")
cat("Sample size:", nrow(df), "\n")
cat("Number of areas:", length(unique(df$area)), "\n")
cat("Number of schools:", length(unique(df$school)), "\n")
cat("Average students per school:", round(nrow(df)/length(unique(df$school)), 2), "\n")
cat("Average schools per area:", round(length(unique(df$school))/length(unique(df$area)), 2), "\n\n")

# Check hierarchical structure
area_school_summary <- df %>%
  group_by(area) %>%
  summarise(
    n_schools = n_distinct(school),
    n_students = n(),
    .groups = 'drop'
  )

school_summary <- df %>%
  group_by(school) %>%
  summarise(
    area = first(area),
    n_students = n(),
    .groups = 'drop'
  )

cat("Areas with number of schools and students:\n")
print(area_school_summary)
cat("\n")

# ===============================================================================
# PREPARE DATA FOR MODELING
# ===============================================================================

# Remove rows with missing values for key variables
model_vars <- c("delinquency", "mgb", "stegb", "ltegb", "vbs", "nrwa", "nrws", 
                "nrwp", "pf", "np", "nle", "sne", "apc", "ats", "atp", "D45", 
                "D47", "bcv", "ica", "delinf", "fses", "age", "gender", "lbc0",
                "school", "area")

# Check for missing values
missing_summary <- df %>%
  select(all_of(model_vars)) %>%
  summarise_all(~sum(is.na(.)))

cat("Missing values by variable:\n")
print(missing_summary)
cat("\n")

# Create analysis dataset with complete cases
df_model <- df %>%
  select(all_of(model_vars)) %>%
  na.omit()

cat("Analysis sample size after removing missing values:", nrow(df_model), "\n")
cat("Percentage of original data retained:", 
    round(100 * nrow(df_model) / nrow(df), 2), "%\n\n")

# Ensure proper factor coding
df_model <- df_model %>%
  mutate(
    gender = factor(gender, levels = c(0, 1), labels = c("Male", "Female")),
    lbc0 = factor(lbc0, levels = c(0, 1, 2), 
                  labels = c("Non-LBC", "Father-Mother-migrant", "Both-migrant")),
    school = factor(school),
    area = factor(area)
  )

# ===============================================================================
# MODEL SPECIFICATION
# ===============================================================================

# Level 1 fixed effects (main effects only - no interactions)
# Base model: mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
#             nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
#             fses + age + gender + lbc0

formula_main <- "delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
                 nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
                 fses + age + gender + lbc0"

# Add random effects for three-level structure
formula_full <- paste0(formula_main, " + (1|area/school)")

cat("=== MODEL FORMULA ===\n")
cat("Fixed effects (main effects only):\n")
cat(formula_main, "\n\n")
cat("Full three-level model:\n") 
cat(formula_full, "\n\n")

# ===============================================================================
# FIT NULL THREE-LEVEL MODEL (BASELINE)
# ===============================================================================

cat("=== FITTING NULL THREE-LEVEL MODEL (BASELINE) ===\n")
cat("This model includes only random effects to establish baseline variance partitioning...\n\n")

# Null model formula (intercept only with random effects)
formula_null <- "delinquency ~ 1 + (1|area/school)"

cat("Null model formula:\n")
cat(formula_null, "\n\n")

# Fit the null three-level negative binomial model
model_null_3level <- glmmTMB(
  formula = as.formula(formula_null),
  data = df_model,
  family = nbinom2(),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

# Null model summary
cat("=== NULL MODEL SUMMARY ===\n")
summary(model_null_3level)
cat("\n")

# ===============================================================================
# NULL MODEL DIAGNOSTICS AND ICC CALCULATION
# ===============================================================================

cat("=== NULL MODEL DIAGNOSTICS ===\n")

# Check convergence
if(model_null_3level$fit$convergence == 0) {
  cat("✓ Null model converged successfully\n")
} else {
  cat("⚠ Warning: Null model convergence issues\n")
}

# Null model fit statistics
cat("\nNull model fit statistics:\n")
cat("AIC:", AIC(model_null_3level), "\n")
cat("BIC:", BIC(model_null_3level), "\n")
cat("Log-likelihood:", logLik(model_null_3level), "\n")

# Random effects from null model
cat("\n=== NULL MODEL RANDOM EFFECTS ===\n")
random_effects_null <- summary(model_null_3level)$varcor
print(random_effects_null)

# Calculate ICC from null model
var_area_null <- as.numeric(random_effects_null$cond$area[1])
var_school_null <- as.numeric(random_effects_null$cond$`school:area`[1])
var_residual_null <- pi^2/3  # For negative binomial, residual variance approximation

total_var_null <- var_area_null + var_school_null + var_residual_null
icc_area_null <- var_area_null / total_var_null
icc_school_null <- var_school_null / total_var_null

cat("\n=== BASELINE INTRACLASS CORRELATIONS (FROM NULL MODEL) ===\n")
cat("Area variance (Level 3):", round(var_area_null, 4), "\n")
cat("School variance (Level 2):", round(var_school_null, 4), "\n")
cat("Residual variance (Level 1):", round(var_residual_null, 4), "\n")
cat("Total variance:", round(total_var_null, 4), "\n\n")

cat("ICC Area (Level 3):", round(icc_area_null, 4), "\n")
cat("ICC School (Level 2):", round(icc_school_null, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area_null + var_school_null) / total_var_null, 4), "\n\n")

cat("Interpretation:\n")
cat("- ICC Area: Proportion of variance between areas\n")
cat("- ICC School: Proportion of variance between schools within areas\n") 
cat("- ICC Combined: Proportion of variance explained by clustering (areas + schools)\n")
cat("- Remaining variance (", round(var_residual_null/total_var_null, 4), ") is at individual level\n\n")

# ===============================================================================
# FIT FULL THREE-LEVEL MODEL WITH PREDICTORS
# ===============================================================================

cat("=== FITTING FULL THREE-LEVEL NEGATIVE BINOMIAL MODEL (WITH PREDICTORS) ===\n")
cat("This may take several minutes...\n\n")

# Fit the three-level negative binomial model using glmmTMB
model_3level_lbc0 <- glmmTMB(
  formula = as.formula(formula_full),
  data = df_model,
  family = nbinom2(),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

# Model summary
cat("=== MODEL SUMMARY ===\n")
summary(model_3level_lbc0)
cat("\n")

# ===============================================================================
# MODEL DIAGNOSTICS
# ===============================================================================

cat("=== MODEL DIAGNOSTICS ===\n")

# Check convergence
if(model_3level_lbc0$fit$convergence == 0) {
  cat("✓ Model converged successfully\n")
} else {
  cat("⚠ Warning: Model convergence issues\n")
}

# Model fit statistics
cat("\nModel fit statistics:\n")
cat("AIC:", AIC(model_3level_lbc0), "\n")
cat("BIC:", BIC(model_3level_lbc0), "\n")
cat("Log-likelihood:", logLik(model_3level_lbc0), "\n")

# Random effects summary
cat("\n=== RANDOM EFFECTS ===\n")
random_effects <- summary(model_3level_lbc0)$varcor
print(random_effects)

# Calculate ICC
var_area <- as.numeric(random_effects$cond$area[1])
var_school <- as.numeric(random_effects$cond$`school:area`[1])
var_residual <- pi^2/3  # For negative binomial, residual variance approximation

total_var <- var_area + var_school + var_residual
icc_area <- var_area / total_var
icc_school <- var_school / total_var

cat("\nIntraclass Correlations:\n")
cat("ICC Area (Level 3):", round(icc_area, 4), "\n")
cat("ICC School (Level 2):", round(icc_school, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area + var_school) / total_var, 4), "\n")

# ===============================================================================
# FIXED EFFECTS RESULTS
# ===============================================================================

cat("\n=== FIXED EFFECTS RESULTS ===\n")

# Extract fixed effects
fixed_effects <- broom.mixed::tidy(model_3level_lbc0, effects = "fixed")
print(fixed_effects)

# Exponentiate coefficients for interpretation (incident rate ratios)
fixed_effects_exp <- fixed_effects %>%
  mutate(
    IRR = exp(estimate),
    IRR_lower = exp(estimate - 1.96 * std.error),
    IRR_upper = exp(estimate + 1.96 * std.error)
  )

cat("\nIncident Rate Ratios (IRR) with 95% CI:\n")
print(fixed_effects_exp %>% 
      select(term, IRR, IRR_lower, IRR_upper, p.value) %>%
      mutate(across(c(IRR, IRR_lower, IRR_upper), ~round(.x, 3))))

# ===============================================================================
# MAIN EFFECTS ANALYSIS
# ===============================================================================

cat("\n=== MAIN EFFECTS ANALYSIS ===\n")

# Extract significant main effects
significant_effects <- fixed_effects_exp %>%
  filter(p.value < 0.05 & term != "(Intercept)") %>%
  arrange(p.value)

cat("Significant main effects (p < 0.05):\n")
if(nrow(significant_effects) > 0) {
  print(significant_effects %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value))
} else {
  cat("No significant main effects found.\n")
}

# ===============================================================================
# MODEL COMPARISONS
# ===============================================================================

cat("\n=== MODEL COMPARISONS ===\n")

# Fit single-level negative binomial model for comparison
model_1level_lbc0 <- glmmTMB(
  formula = as.formula(formula_main),
  data = df_model,
  family = nbinom2()
)

cat("Model comparison (AIC):\n")
cat("Null three-level model AIC:", AIC(model_null_3level), "\n")
cat("Single-level model AIC:", AIC(model_1level_lbc0), "\n")
cat("Full three-level model AIC:", AIC(model_3level_lbc0), "\n\n")

cat("AIC improvements:\n")
cat("Three-level vs Single-level:", AIC(model_1level_lbc0) - AIC(model_3level_lbc0), "\n")
cat("Full vs Null three-level:", AIC(model_null_3level) - AIC(model_3level_lbc0), "\n\n")

# Likelihood ratio tests
cat("Likelihood ratio tests:\n")
lrt_null_vs_full <- anova(model_null_3level, model_3level_lbc0)
cat("Null vs Full three-level model:\n")
print(lrt_null_vs_full)

lrt_single_vs_full <- anova(model_1level_lbc0, model_3level_lbc0)
cat("\nSingle-level vs Full three-level model:\n")
print(lrt_single_vs_full)

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save model objects
saveRDS(model_null_3level, "model_null_3level_nb.rds")
saveRDS(model_3level_lbc0, "model_3level_nb_lbc0.rds")

# Save fixed effects results
write_csv(fixed_effects_exp, "fixed_effects_results_lbc0.csv")

# Save comprehensive model summary
sink("model_summary_lbc0.txt")
cat("Three-Level Hierarchical Negative Binomial Regression Results (Simplified - LBC0)\n")
cat("=============================================================================\n\n")
cat("Sample size:", nrow(df_model), "\n")
cat("Number of areas:", length(unique(df_model$area)), "\n")
cat("Number of schools:", length(unique(df_model$school)), "\n\n")

cat("=== NULL MODEL RESULTS ===\n")
cat("Baseline Intraclass Correlations (from null model):\n")
cat("ICC Area (Level 3):", round(icc_area_null, 4), "\n")
cat("ICC School (Level 2):", round(icc_school_null, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area_null + var_school_null) / total_var_null, 4), "\n\n")

cat("=== FULL MODEL RESULTS ===\n")
summary(model_3level_lbc0)
cat("\n\nFull Model Intraclass Correlations:\n")
cat("ICC Area (Level 3):", round(icc_area, 4), "\n")
cat("ICC School (Level 2):", round(icc_school, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area + var_school) / total_var, 4), "\n")

cat("\n=== MODEL COMPARISONS ===\n")
cat("AIC Comparisons:\n")
cat("Null three-level model:", AIC(model_null_3level), "\n")
cat("Full three-level model:", AIC(model_3level_lbc0), "\n")
cat("AIC improvement (predictors):", AIC(model_null_3level) - AIC(model_3level_lbc0), "\n")
sink()

# Save null model ICC results separately
null_icc_results <- data.frame(
  Component = c("Area (Level 3)", "School (Level 2)", "Combined (Levels 2+3)", "Individual (Level 1)"),
  Variance = c(var_area_null, var_school_null, var_area_null + var_school_null, var_residual_null),
  ICC = c(icc_area_null, icc_school_null, (var_area_null + var_school_null) / total_var_null, var_residual_null / total_var_null)
)
write_csv(null_icc_results, "null_model_icc_results.csv")

cat("Results saved:\n")
cat("- Null model object: model_null_3level_nb.rds\n")
cat("- Full model object: model_3level_nb_lbc0.rds\n")
cat("- Fixed effects: fixed_effects_results_lbc0.csv\n")
cat("- Null model ICC: null_model_icc_results.csv\n")
cat("- Comprehensive summary: model_summary_lbc0.txt\n")

cat("\n=== SIMPLIFIED ANALYSIS COMPLETE (LBC0) ===\n") 