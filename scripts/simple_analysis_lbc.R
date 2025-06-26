# ===============================================================================
# Three-Level Hierarchical Negative Binomial Regression Analysis (Simplified + Complex)
# ICP-H Project - Using lbc instead of lbc0
# 
# Level 1: Individual (main effects + interactions with lbc)
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
                "D47", "bcv", "ica", "delinf", "fses", "age", "gender", "lbc",
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
    lbc = factor(lbc, levels = c(0, 1), 
                 labels = c("Non-LBC", "LBC")),
    school = factor(school),
    area = factor(area)
  )

# ===============================================================================
# MODEL SPECIFICATION
# ===============================================================================

# Level 1 fixed effects - SIMPLE MODEL (main effects only)
formula_main_simple <- "delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
                        nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
                        fses + age + gender + lbc"

# Level 1 fixed effects - COMPLEX MODEL (main effects + lbc interactions)
formula_main_complex <- "delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
                         nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
                         fses + age + gender + lbc + 
                         lbc:mgb + lbc:stegb + lbc:ltegb + lbc:vbs + lbc:nrwa + lbc:nrwp + 
                         lbc:pf + lbc:np + lbc:nle"

# Add random effects for three-level structure
formula_simple_full <- paste0(formula_main_simple, " + (1|area/school)")
formula_complex_full <- paste0(formula_main_complex, " + (1|area/school)")

cat("=== MODEL FORMULAS ===\n")
cat("Simple model (main effects only):\n")
cat(formula_main_simple, "\n\n")
cat("Complex model (main effects + lbc interactions):\n")
cat(formula_main_complex, "\n\n")
cat("Simple three-level model:\n") 
cat(formula_simple_full, "\n\n")
cat("Complex three-level model:\n") 
cat(formula_complex_full, "\n\n")

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
# FIT SIMPLE THREE-LEVEL MODEL WITH MAIN EFFECTS ONLY
# ===============================================================================

cat("=== FITTING SIMPLE THREE-LEVEL NEGATIVE BINOMIAL MODEL (MAIN EFFECTS ONLY) ===\n")
cat("This may take several minutes...\n\n")

# Fit the simple three-level negative binomial model using glmmTMB
model_3level_simple <- glmmTMB(
  formula = as.formula(formula_simple_full),
  data = df_model,
  family = nbinom2(),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

# Model summary
cat("=== SIMPLE MODEL SUMMARY ===\n")
summary(model_3level_simple)
cat("\n")

# ===============================================================================
# FIT COMPLEX THREE-LEVEL MODEL WITH INTERACTIONS
# ===============================================================================

cat("=== FITTING COMPLEX THREE-LEVEL NEGATIVE BINOMIAL MODEL (WITH LBC INTERACTIONS) ===\n")
cat("This may take several minutes (longer than simple model)...\n\n")

# Fit the complex three-level negative binomial model using glmmTMB
model_3level_complex <- glmmTMB(
  formula = as.formula(formula_complex_full),
  data = df_model,
  family = nbinom2(),
  control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
)

# Model summary
cat("=== COMPLEX MODEL SUMMARY ===\n")
summary(model_3level_complex)
cat("\n")

# ===============================================================================
# MODEL DIAGNOSTICS
# ===============================================================================

cat("=== MODEL DIAGNOSTICS ===\n")

# Check convergence for both models
if(model_3level_simple$fit$convergence == 0) {
  cat("✓ Simple model converged successfully\n")
} else {
  cat("⚠ Warning: Simple model convergence issues\n")
}

if(model_3level_complex$fit$convergence == 0) {
  cat("✓ Complex model converged successfully\n")
} else {
  cat("⚠ Warning: Complex model convergence issues\n")
}

# Model fit statistics comparison
cat("\nModel fit statistics:\n")
cat("Simple model AIC:", AIC(model_3level_simple), "\n")
cat("Complex model AIC:", AIC(model_3level_complex), "\n")
cat("AIC improvement (complex vs simple):", AIC(model_3level_simple) - AIC(model_3level_complex), "\n\n")

cat("Simple model BIC:", BIC(model_3level_simple), "\n")
cat("Complex model BIC:", BIC(model_3level_complex), "\n")
cat("BIC improvement (complex vs simple):", BIC(model_3level_simple) - BIC(model_3level_complex), "\n\n")

# Random effects summary for complex model
cat("\n=== COMPLEX MODEL RANDOM EFFECTS ===\n")
random_effects_complex <- summary(model_3level_complex)$varcor
print(random_effects_complex)

# Calculate ICC for complex model
var_area_complex <- as.numeric(random_effects_complex$cond$area[1])
var_school_complex <- as.numeric(random_effects_complex$cond$`school:area`[1])
var_residual_complex <- pi^2/3  # For negative binomial, residual variance approximation

total_var_complex <- var_area_complex + var_school_complex + var_residual_complex
icc_area_complex <- var_area_complex / total_var_complex
icc_school_complex <- var_school_complex / total_var_complex

cat("\nComplex Model Intraclass Correlations:\n")
cat("ICC Area (Level 3):", round(icc_area_complex, 4), "\n")
cat("ICC School (Level 2):", round(icc_school_complex, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area_complex + var_school_complex) / total_var_complex, 4), "\n")

# ===============================================================================
# FIXED EFFECTS RESULTS
# ===============================================================================

cat("\n=== FIXED EFFECTS RESULTS ===\n")

# Extract fixed effects for both models
fixed_effects_simple <- broom.mixed::tidy(model_3level_simple, effects = "fixed")
fixed_effects_complex <- broom.mixed::tidy(model_3level_complex, effects = "fixed")

cat("Simple Model Fixed Effects:\n")
print(fixed_effects_simple)
cat("\n")

cat("Complex Model Fixed Effects:\n")
print(fixed_effects_complex)
cat("\n")

# Exponentiate coefficients for interpretation (incident rate ratios)
fixed_effects_simple_exp <- fixed_effects_simple %>%
  mutate(
    IRR = exp(estimate),
    IRR_lower = exp(estimate - 1.96 * std.error),
    IRR_upper = exp(estimate + 1.96 * std.error)
  )

fixed_effects_complex_exp <- fixed_effects_complex %>%
  mutate(
    IRR = exp(estimate),
    IRR_lower = exp(estimate - 1.96 * std.error),
    IRR_upper = exp(estimate + 1.96 * std.error)
  )

cat("\nSimple Model - Incident Rate Ratios (IRR) with 95% CI:\n")
print(fixed_effects_simple_exp %>% 
      select(term, IRR, IRR_lower, IRR_upper, p.value) %>%
      mutate(across(c(IRR, IRR_lower, IRR_upper), ~round(.x, 3))))

cat("\nComplex Model - Incident Rate Ratios (IRR) with 95% CI:\n")
print(fixed_effects_complex_exp %>% 
      select(term, IRR, IRR_lower, IRR_upper, p.value) %>%
      mutate(across(c(IRR, IRR_lower, IRR_upper), ~round(.x, 3))))

# ===============================================================================
# INTERACTION EFFECTS ANALYSIS
# ===============================================================================

cat("\n=== INTERACTION EFFECTS ANALYSIS ===\n")

# Extract significant interaction effects from complex model
interaction_effects <- fixed_effects_complex_exp %>%
  filter(grepl("lbc", term) & term != "lbcLBC") %>%  # Get interaction terms
  arrange(p.value)

cat("LBC Interaction Effects (from Complex Model):\n")
if(nrow(interaction_effects) > 0) {
  print(interaction_effects %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value) %>%
        mutate(across(c(IRR, IRR_lower, IRR_upper), ~round(.x, 3))))
} else {
  cat("No interaction effects found.\n")
}

# Identify significant interactions
significant_interactions <- interaction_effects %>%
  filter(p.value < 0.05)

cat("\nSignificant LBC Interaction Effects (p < 0.05):\n")
if(nrow(significant_interactions) > 0) {
  print(significant_interactions %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value))
} else {
  cat("No significant interaction effects found.\n")
}

# ===============================================================================
# MAIN EFFECTS ANALYSIS
# ===============================================================================

cat("\n=== MAIN EFFECTS ANALYSIS ===\n")

# Extract significant main effects from complex model
significant_main_effects <- fixed_effects_complex_exp %>%
  filter(p.value < 0.05 & term != "(Intercept)" & !grepl(":", term)) %>%
  arrange(p.value)

cat("Significant main effects from Complex Model (p < 0.05):\n")
if(nrow(significant_main_effects) > 0) {
  print(significant_main_effects %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value))
} else {
  cat("No significant main effects found.\n")
}

# ===============================================================================
# MODEL COMPARISONS
# ===============================================================================

cat("\n=== MODEL COMPARISONS ===\n")

# Fit single-level models for comparison
model_1level_simple <- glmmTMB(
  formula = as.formula(formula_main_simple),
  data = df_model,
  family = nbinom2()
)

model_1level_complex <- glmmTMB(
  formula = as.formula(formula_main_complex),
  data = df_model,
  family = nbinom2()
)

cat("Model comparison (AIC):\n")
cat("Null three-level model AIC:", AIC(model_null_3level), "\n")
cat("Single-level simple model AIC:", AIC(model_1level_simple), "\n")
cat("Single-level complex model AIC:", AIC(model_1level_complex), "\n")
cat("Three-level simple model AIC:", AIC(model_3level_simple), "\n")
cat("Three-level complex model AIC:", AIC(model_3level_complex), "\n\n")

cat("AIC improvements:\n")
cat("Three-level simple vs Single-level simple:", AIC(model_1level_simple) - AIC(model_3level_simple), "\n")
cat("Three-level complex vs Single-level complex:", AIC(model_1level_complex) - AIC(model_3level_complex), "\n")
cat("Complex vs Simple (three-level):", AIC(model_3level_simple) - AIC(model_3level_complex), "\n")
cat("Complex vs Null (three-level):", AIC(model_null_3level) - AIC(model_3level_complex), "\n\n")

# Likelihood ratio tests
cat("Likelihood ratio tests:\n")
lrt_null_vs_simple <- anova(model_null_3level, model_3level_simple)
cat("Null vs Simple three-level model:\n")
print(lrt_null_vs_simple)

lrt_simple_vs_complex <- anova(model_3level_simple, model_3level_complex)
cat("\nSimple vs Complex three-level model:\n")
print(lrt_simple_vs_complex)

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

cat("\n=== SAVING RESULTS ===\n")

# Save model objects
saveRDS(model_null_3level, "model_null_3level_nb_lbc.rds")
saveRDS(model_3level_simple, "model_3level_simple_nb_lbc.rds")
saveRDS(model_3level_complex, "model_3level_complex_nb_lbc.rds")

# Save fixed effects results
write_csv(fixed_effects_simple_exp, "fixed_effects_simple_lbc.csv")
write_csv(fixed_effects_complex_exp, "fixed_effects_complex_lbc.csv")
write_csv(interaction_effects, "interaction_effects_lbc.csv")

# Save comprehensive model summary
sink("model_summary_complex_lbc.txt")
cat("Three-Level Hierarchical Negative Binomial Regression Results (Simple + Complex - LBC)\n")
cat("========================================================================================\n\n")
cat("Sample size:", nrow(df_model), "\n")
cat("Number of areas:", length(unique(df_model$area)), "\n")
cat("Number of schools:", length(unique(df_model$school)), "\n\n")

cat("=== NULL MODEL RESULTS ===\n")
cat("Baseline Intraclass Correlations (from null model):\n")
cat("ICC Area (Level 3):", round(icc_area_null, 4), "\n")
cat("ICC School (Level 2):", round(icc_school_null, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area_null + var_school_null) / total_var_null, 4), "\n\n")

cat("=== SIMPLE MODEL RESULTS ===\n")
summary(model_3level_simple)

cat("\n=== COMPLEX MODEL RESULTS ===\n")
summary(model_3level_complex)

cat("\n\nComplex Model Intraclass Correlations:\n")
cat("ICC Area (Level 3):", round(icc_area_complex, 4), "\n")
cat("ICC School (Level 2):", round(icc_school_complex, 4), "\n")
cat("ICC Combined (Levels 2+3):", round((var_area_complex + var_school_complex) / total_var_complex, 4), "\n")

cat("\n=== MODEL COMPARISONS ===\n")
cat("AIC Comparisons:\n")
cat("Null three-level model:", AIC(model_null_3level), "\n")
cat("Simple three-level model:", AIC(model_3level_simple), "\n")
cat("Complex three-level model:", AIC(model_3level_complex), "\n")
cat("AIC improvement (simple vs null):", AIC(model_null_3level) - AIC(model_3level_simple), "\n")
cat("AIC improvement (complex vs simple):", AIC(model_3level_simple) - AIC(model_3level_complex), "\n")

cat("\n=== SIGNIFICANT INTERACTION EFFECTS ===\n")
if(nrow(significant_interactions) > 0) {
  print(significant_interactions %>% 
        select(term, IRR, IRR_lower, IRR_upper, p.value))
} else {
  cat("No significant interaction effects found.\n")
}
sink()

# Save interaction effects analysis
interaction_summary <- data.frame(
  Model = "Complex",
  N_Interactions = nrow(interaction_effects),
  N_Significant_Interactions = nrow(significant_interactions),
  AIC_Simple = AIC(model_3level_simple),
  AIC_Complex = AIC(model_3level_complex),
  AIC_Improvement = AIC(model_3level_simple) - AIC(model_3level_complex)
)
write_csv(interaction_summary, "interaction_summary_lbc.csv")

cat("Results saved:\n")
cat("- Null model object: model_null_3level_nb_lbc.rds\n")
cat("- Simple model object: model_3level_simple_nb_lbc.rds\n")
cat("- Complex model object: model_3level_complex_nb_lbc.rds\n")
cat("- Simple model fixed effects: fixed_effects_simple_lbc.csv\n")
cat("- Complex model fixed effects: fixed_effects_complex_lbc.csv\n")
cat("- Interaction effects: interaction_effects_lbc.csv\n")
cat("- Interaction summary: interaction_summary_lbc.csv\n")
cat("- Comprehensive summary: model_summary_complex_lbc.txt\n")

cat("\n=== ANALYSIS COMPLETE (SIMPLE + COMPLEX MODELS WITH LBC) ===\n")
