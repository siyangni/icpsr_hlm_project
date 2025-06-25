# ===============================================================================
# Three-Level Hierarchical Negative Binomial Regression Analysis (Simplified)
# ICP-H Project
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
                "D47", "bcv", "ica", "delinf", "fses", "age", "gender", "lbc1",
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
    lbc1 = factor(lbc1, levels = c(0, 1, 2, 3), 
                  labels = c("Non-LBC", "Father-migrant", "Mother-migrant", "Both-migrant")),
    school = factor(school),
    area = factor(area)
  )

# ===============================================================================
# MODEL SPECIFICATION
# ===============================================================================

# Level 1 fixed effects (main effects only - no interactions)
# Base model: mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
#             nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
#             fses + age + gender + lbc1

formula_main <- "delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
                 nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
                 fses + age + gender + lbc1"

# Add random effects for three-level structure
formula_full <- paste0(formula_main, " + (1|area/school)")

cat("=== MODEL FORMULA ===\n")
cat("Fixed effects (main effects only):\n")
cat(formula_main, "\n\n")
cat("Full three-level model:\n") 
cat(formula_full, "\n\n")

# ===============================================================================
# FIT HIERARCHICAL NEGATIVE BINOMIAL MODEL
# ===============================================================================

cat("=== FITTING THREE-LEVEL NEGATIVE BINOMIAL MODEL (SIMPLIFIED) ===\n")
cat("This may take several minutes...\n\n")

# Fit the three-level negative binomial model using glmmTMB
model_3level_simple <- glmmTMB(
  formula = as.formula(formula_full),
  data = df_model,
  family = nbinom2(),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

# Model summary
cat("=== MODEL SUMMARY ===\n")
summary(model_3level_simple)
cat("\n")

# ===============================================================================
# MODEL DIAGNOSTICS
# ===============================================================================

cat("=== MODEL DIAGNOSTICS ===\n")

# Check convergence
if(model_3level_simple$fit$convergence == 0) {
  cat("✓ Model converged successfully\n")
} else {
  cat("⚠ Warning: Model convergence issues\n")
}

# Model fit statistics
cat("\nModel fit statistics:\n")
cat("AIC:", AIC(model_3level_simple), "\n")
cat("BIC:", BIC(model_3level_simple), "\n")
cat("Log-likelihood:", logLik(model_3level_simple), "\n")

# Random effects summary
cat("\n=== RANDOM EFFECTS ===\n")
random_effects <- summary(model_3level_simple)$varcor
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
fixed_effects <- broom.mixed::tidy(model_3level_simple, effects = "fixed")
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
# COMPARISON WITH SINGLE-LEVEL MODEL
# ===============================================================================

cat("\n=== COMPARISON WITH SINGLE-LEVEL MODEL ===\n")

# Fit single-level negative binomial model for comparison
model_1level_simple <- glmmTMB(
  formula = as.formula(formula_main),
  data = df_model,
  family = nbinom2()
)

cat("Single-level model AIC:", AIC(model_1level_simple), "\n")
cat("Three-level model AIC:", AIC(model_3level_simple), "\n")
cat("AIC improvement:", AIC(model_1level_simple) - AIC(model_3level_simple), "\n")

# Likelihood ratio test
lrt_result <- anova(model_1level_simple, model_3level_simple)
print(lrt_result)

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save model object
saveRDS(model_3level_simple, "model_3level_nb_simple.rds")

# Save fixed effects results
write_csv(fixed_effects_exp, "fixed_effects_results_simple.csv")

# Save model summary
sink("model_summary_simple.txt")
cat("Three-Level Hierarchical Negative Binomial Regression Results (Simplified)\n")
cat("=======================================================================\n\n")
cat("Sample size:", nrow(df_model), "\n")
cat("Number of areas:", length(unique(df_model$area)), "\n")
cat("Number of schools:", length(unique(df_model$school)), "\n\n")
summary(model_3level_simple)
cat("\n\nIntraclass Correlations:\n")
cat("ICC Area (Level 3):", round(icc_area, 4), "\n")
cat("ICC School (Level 2):", round(icc_school, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area + var_school) / total_var, 4), "\n")
sink()

cat("Results saved:\n")
cat("- Model object: model_3level_nb_simple.rds\n")
cat("- Fixed effects: fixed_effects_results_simple.csv\n") 
cat("- Full summary: model_summary_simple.txt\n")

cat("\n=== SIMPLIFIED ANALYSIS COMPLETE ===\n") 