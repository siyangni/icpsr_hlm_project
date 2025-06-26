# ===============================================================================
# Three-Level Hierarchical Negative Binomial Regression Analysis
# ICP-H Project
# 
# Level 1: Individual (with interaction effects)
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

# Level 1 fixed effects (from ch4_regression.R lines 103-107 plus interactions)
# Base model: mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
#             nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
#             fses + age + gender + lbc1
# 
# Interaction effects: lbc1 * (mgb, stegb, ltegb, vbs, nrwa, nrwp, pf, np, nle)

formula_base <- "delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
                 nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
                 fses + age + gender + lbc1"

# Add interaction terms
formula_interactions <- paste0(formula_base, " + ",
                              "lbc1:mgb + lbc1:stegb + lbc1:ltegb + lbc1:vbs + ",
                              "lbc1:nrwa + lbc1:nrwp + lbc1:pf + lbc1:np + lbc1:nle")

# Add random effects for three-level structure
formula_full <- paste0(formula_interactions, " + (1|area/school)")

cat("=== MODEL FORMULA ===\n")
cat("Fixed effects with interactions:\n")
cat(formula_interactions, "\n\n")
cat("Full three-level model:\n") 
cat(formula_full, "\n\n")

# ===============================================================================
# FIT HIERARCHICAL NEGATIVE BINOMIAL MODEL
# ===============================================================================

cat("=== FITTING THREE-LEVEL NEGATIVE BINOMIAL MODEL ===\n")
cat("This may take several minutes...\n\n")

# Fit the three-level negative binomial model using glmmTMB
model_3level <- glmmTMB(
  formula = as.formula(formula_full),
  data = df_model,
  family = nbinom2(),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

# Model summary
cat("=== MODEL SUMMARY ===\n")
summary(model_3level)
cat("\n")

# ===============================================================================
# MODEL DIAGNOSTICS
# ===============================================================================

cat("=== MODEL DIAGNOSTICS ===\n")

# Check convergence
if(model_3level$fit$convergence == 0) {
  cat("✓ Model converged successfully\n")
} else {
  cat("⚠ Warning: Model convergence issues\n")
}

# Model fit statistics
cat("\nModel fit statistics:\n")
cat("AIC:", AIC(model_3level), "\n")
cat("BIC:", BIC(model_3level), "\n")
cat("Log-likelihood:", logLik(model_3level), "\n")

# Random effects summary
cat("\n=== RANDOM EFFECTS ===\n")
random_effects <- summary(model_3level)$varcor
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
fixed_effects <- broom.mixed::tidy(model_3level, effects = "fixed")
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
# INTERACTION EFFECTS ANALYSIS
# ===============================================================================

cat("\n=== INTERACTION EFFECTS ANALYSIS ===\n")

# Extract interaction terms
interaction_terms <- fixed_effects_exp %>%
  filter(grepl("lbc1.*:", term)) %>%
  arrange(p.value)

cat("Significant interaction effects (p < 0.05):\n")
significant_interactions <- interaction_terms %>%
  filter(p.value < 0.05)

if(nrow(significant_interactions) > 0) {
  print(significant_interactions %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value))
} else {
  cat("No significant interaction effects found.\n")
}

cat("\nAll interaction effects:\n")
print(interaction_terms %>% 
      select(term, IRR, IRR_lower, IRR_upper, p.value))

# ===============================================================================
# COMPARISON WITH SINGLE-LEVEL MODEL
# ===============================================================================

cat("\n=== COMPARISON WITH SINGLE-LEVEL MODEL ===\n")

# Fit single-level negative binomial model for comparison
model_1level <- glmmTMB(
  formula = as.formula(formula_interactions),
  data = df_model,
  family = nbinom2()
)

cat("Single-level model AIC:", AIC(model_1level), "\n")
cat("Three-level model AIC:", AIC(model_3level), "\n")
cat("AIC improvement:", AIC(model_1level) - AIC(model_3level), "\n")

# Likelihood ratio test
lrt_result <- anova(model_1level, model_3level)
print(lrt_result)

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save model object
saveRDS(model_3level, "model_3level_nb.rds")

# Save fixed effects results
write_csv(fixed_effects_exp, "fixed_effects_results.csv")

# Save model summary
sink("model_summary.txt")
cat("Three-Level Hierarchical Negative Binomial Regression Results\n")
cat("============================================================\n\n")
cat("Sample size:", nrow(df_model), "\n")
cat("Number of areas:", length(unique(df_model$area)), "\n")
cat("Number of schools:", length(unique(df_model$school)), "\n\n")
summary(model_3level)
cat("\n\nIntraclass Correlations:\n")
cat("ICC Area (Level 3):", round(icc_area, 4), "\n")
cat("ICC School (Level 2):", round(icc_school, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area + var_school) / total_var, 4), "\n")
sink()

cat("Results saved:\n")
cat("- Model object: model_3level_nb.rds\n")
cat("- Fixed effects: fixed_effects_results.csv\n") 
cat("- Full summary: model_summary.txt\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
