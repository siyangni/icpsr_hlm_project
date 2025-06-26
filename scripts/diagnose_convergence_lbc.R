# ===============================================================================
# Diagnostic Script for Convergence Issues in simple_analysis_lbc.R
# ===============================================================================

library(pacman)
p_load(dplyr, readr, glmmTMB, broom.mixed)

# Load the data and models
df <- read_csv("data/processed.csv")
model_null_3level <- readRDS("model_null_3level_nb_lbc.rds")
model_3level_lbc <- readRDS("model_3level_nb_lbc.rds")

cat("=== CONVERGENCE DIAGNOSTICS ===\n")

# 1. Check convergence status
cat("1. CONVERGENCE STATUS:\n")
cat("Null model convergence code:", model_null_3level$fit$convergence, "\n")
cat("Full model convergence code:", model_3level_lbc$fit$convergence, "\n")
cat("(0 = successful convergence)\n\n")

# 2. Check gradient and Hessian
cat("2. GRADIENT AND HESSIAN DIAGNOSTICS:\n")
cat("Full model gradient norm:", max(abs(model_3level_lbc$fit$par)), "\n")

# Check if Hessian is positive definite
hessian_eigen <- try(eigen(model_3level_lbc$fit$hessian, only.values = TRUE)$values)
if(!inherits(hessian_eigen, "try-error")) {
  cat("Minimum Hessian eigenvalue:", min(hessian_eigen), "\n")
  cat("Hessian is positive definite:", all(hessian_eigen > 0), "\n")
} else {
  cat("Could not compute Hessian eigenvalues\n")
}
cat("\n")

# 3. Check random effects variances
cat("3. RANDOM EFFECTS VARIANCE DIAGNOSTICS:\n")
null_varcor <- VarCorr(model_null_3level)
full_varcor <- VarCorr(model_3level_lbc)

cat("Null model variances:\n")
print(null_varcor)
cat("\nFull model variances:\n")
print(full_varcor)

# Extract variance components
null_area_var <- as.numeric(attr(null_varcor$cond$area, "stddev"))^2
null_school_var <- as.numeric(attr(null_varcor$cond$`school:area`, "stddev"))^2

full_area_var <- as.numeric(attr(full_varcor$cond$area, "stddev"))^2
full_school_var <- as.numeric(attr(full_varcor$cond$`school:area`, "stddev"))^2

cat("\nVariance components:\n")
cat("Null model - Area variance:", null_area_var, "School variance:", null_school_var, "\n")
cat("Full model - Area variance:", full_area_var, "School variance:", full_school_var, "\n")

# Check if school variance is near zero (boundary issue)
if(full_school_var < 0.01) {
  cat("⚠ WARNING: School-level variance is very small (", full_school_var, ")\n")
  cat("This suggests the three-level structure may not be necessary.\n")
}
cat("\n")

# 4. Data structure diagnostics
cat("4. DATA STRUCTURE DIAGNOSTICS:\n")
model_vars <- c("delinquency", "mgb", "stegb", "ltegb", "vbs", "nrwa", "nrws", 
                "nrwp", "pf", "np", "nle", "sne", "apc", "ats", "atp", "D45", 
                "D47", "bcv", "ica", "delinf", "fses", "age", "gender", "lbc",
                "school", "area")

df_model <- df %>%
  select(all_of(model_vars)) %>%
  na.omit() %>%
  mutate(
    gender = factor(gender, levels = c(0, 1), labels = c("Male", "Female")),
    lbc = factor(lbc, levels = c(0, 1), labels = c("Non-LBC", "LBC")),
    school = factor(school),
    area = factor(area)
  )

# Check hierarchical structure
school_counts <- df_model %>%
  group_by(area, school) %>%
  summarise(n = n(), .groups = 'drop')

cat("Students per school (summary):\n")
print(summary(school_counts$n))

cat("\nSchools with very few students (< 20):\n")
small_schools <- school_counts %>% filter(n < 20)
print(small_schools)

# Check for schools with zero variance in outcome
school_var_check <- df_model %>%
  group_by(school) %>%
  summarise(
    n = n(),
    delinq_var = var(delinquency, na.rm = TRUE),
    delinq_mean = mean(delinquency, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\nSchools with zero variance in delinquency:\n")
zero_var_schools <- school_var_check %>% filter(delinq_var == 0 | is.na(delinq_var))
print(zero_var_schools)

cat("\n=== POTENTIAL SOLUTIONS ===\n")
cat("1. Try a two-level model (remove school level) if school variance is negligible\n")
cat("2. Use different optimizer (e.g., 'bobyqa', 'Nelder_Mead')\n")
cat("3. Scale predictors if they have very different magnitudes\n")
cat("4. Consider removing schools with very few observations\n")
cat("5. Use simpler random effects structure\n\n")

# 5. Try alternative optimizers
cat("5. TESTING ALTERNATIVE OPTIMIZERS:\n")

# Test with different optimizer
cat("Trying bobyqa optimizer...\n")
model_bobyqa <- try({
  glmmTMB(
    formula = delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
              nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
              fses + age + gender + lbc + (1|area/school),
    data = df_model,
    family = nbinom2(),
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "L-BFGS-B"))
  )
})

if(!inherits(model_bobyqa, "try-error")) {
  cat("✓ Alternative optimizer succeeded\n")
  cat("Convergence code:", model_bobyqa$fit$convergence, "\n")
  saveRDS(model_bobyqa, "model_3level_nb_lbc_alt_optimizer.rds")
} else {
  cat("✗ Alternative optimizer also failed\n")
}

# 6. Test two-level model
cat("\nTesting two-level model (area only)...\n")
model_2level <- try({
  glmmTMB(
    formula = delinquency ~ mgb + stegb + ltegb + vbs + nrwa + nrws + nrwp + pf + np + 
              nle + sne + apc + ats + atp + D45 + D47 + bcv + ica + delinf + 
              fses + age + gender + lbc + (1|area),
    data = df_model,
    family = nbinom2()
  )
})

if(!inherits(model_2level, "try-error")) {
  cat("✓ Two-level model converged successfully\n")
  cat("Convergence code:", model_2level$fit$convergence, "\n")
  cat("AIC:", AIC(model_2level), "\n")
  saveRDS(model_2level, "model_2level_nb_lbc.rds")
  
  # Compare with three-level
  cat("AIC comparison - Two-level:", AIC(model_2level), "vs Three-level:", AIC(model_3level_lbc), "\n")
} else {
  cat("✗ Two-level model also failed\n")
}

cat("\n=== COMPLETE MODEL SUMMARY (FIXED) ===\n")

# Generate complete summary with all details
sink("model_summary_lbc_complete.txt")
cat("Three-Level Hierarchical Negative Binomial Regression Results (LBC) - COMPLETE\n")
cat("================================================================================\n\n")

cat("Sample size:", nrow(df_model), "\n")
cat("Number of areas:", length(unique(df_model$area)), "\n")
cat("Number of schools:", length(unique(df_model$school)), "\n\n")

cat("=== CONVERGENCE STATUS ===\n")
cat("Null model convergence:", model_null_3level$fit$convergence, "\n")
cat("Full model convergence:", model_3level_lbc$fit$convergence, "\n")
if(model_3level_lbc$fit$convergence != 0) {
  cat("⚠ WARNING: Convergence issues detected\n")
}
cat("\n")

cat("=== NULL MODEL RESULTS ===\n")
cat("Baseline Intraclass Correlations:\n")
cat("ICC Area (Level 3):", round(null_area_var / (null_area_var + null_school_var + pi^2/3), 4), "\n")
cat("ICC School (Level 2):", round(null_school_var / (null_area_var + null_school_var + pi^2/3), 4), "\n")
cat("ICC Combined:", round((null_area_var + null_school_var) / (null_area_var + null_school_var + pi^2/3), 4), "\n\n")

cat("=== FULL MODEL RESULTS ===\n")
print(summary(model_3level_lbc))

cat("\n=== FIXED EFFECTS (DETAILED) ===\n")
fixed_effects <- broom.mixed::tidy(model_3level_lbc, effects = "fixed")
fixed_effects_exp <- fixed_effects %>%
  mutate(
    IRR = exp(estimate),
    IRR_lower = exp(estimate - 1.96 * std.error),
    IRR_upper = exp(estimate + 1.96 * std.error)
  )
print(fixed_effects_exp)

cat("\n=== MODEL FIT STATISTICS ===\n")
cat("AIC:", AIC(model_3level_lbc), "\n")
cat("BIC:", BIC(model_3level_lbc), "\n")
cat("Log-likelihood:", logLik(model_3level_lbc), "\n")

sink()

cat("Complete summary saved to: model_summary_lbc_complete.txt\n") 