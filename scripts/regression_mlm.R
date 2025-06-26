# Multi-Level Negative Binomial Regression
# Replicating ch4_regression_new.R (lines 60-115) with hierarchical structure
# Level 1: Individual
# Level 2: School  
# Level 3: Area

library(pacman)
p_load(tidyverse, Hmisc, MASS, sandwich, lmtest, coxphw, stargazer, 
       glmmTMB, broom.mixed, sjPlot, performance)

# Functions for robust standard errors and output formatting
res <- function(m){
  cov.m <- vcovHC(m, type = "HC0")
  std.err <- sqrt(diag(cov.m))
  q.val <- qnorm(0.975)
  r.est <- cbind(
    Estimate = coef(m)
    , "Robust SE" = std.err
    , z = (coef(m)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(m)/std.err), lower.tail = FALSE)
    , LL = coef(m) - q.val  * std.err
    , UL = coef(m) + q.val  * std.err)
  print(r.est)}

# Function to extract IRR with confidence intervals from glmmTMB models
get_irr_ci <- function(model) {
  # Extract fixed effects coefficients
  coeff <- fixef(model)$cond
  
  # Get confidence intervals using broom.mixed which handles glmmTMB better
  ci_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  
  # Create IRR table
  irr_table <- data.frame(
    IRR = exp(coeff),
    Lower_CI = exp(ci_df$conf.low),
    Upper_CI = exp(ci_df$conf.high),
    P_value = ci_df$p.value
  )
  
  rownames(irr_table) <- names(coeff)
  return(irr_table)
}

# Function to extract ICC from glmmTMB models
get_icc <- function(model) {
  var_comps <- VarCorr(model)$cond
  var_area <- as.numeric(var_comps$area[1])
  var_school <- as.numeric(var_comps$`school:area`[1])
  
  # For negative binomial, residual variance approximated as pi^2/3
  var_residual <- pi^2/3
  
  icc_area <- var_area / (var_area + var_school + var_residual)
  icc_school <- var_school / (var_area + var_school + var_residual)
  icc_total <- (var_area + var_school) / (var_area + var_school + var_residual)
  
  list(
    icc_area = icc_area,
    icc_school = icc_school, 
    icc_total = icc_total,
    var_area = var_area,
    var_school = var_school
  )
}

# Load the data
load("data/imp_aggregate.RData")
df <- dat$imputations$imp1 # Take the first set of the five

# Preparing the operation data frame for regression
rdf <- df %>%
  dplyr::select(delinquency, mgb, stegb, ltegb, vbs, nrwa,
         nrws, nrwp, pf, np, nle, sne, apc,
         ats, atp, D45, D47, bcv, ica, delinf, fses,
         age, gender, school, lbc, lbc1, area) %>%
  na.omit() %>%
  mutate(
    gender = factor(gender),
    lbc = factor(lbc),
    lbc1 = factor(lbc1),
    school = factor(school),
    area = factor(area)
  )

# Check hierarchical structure
cat("=== DATA OVERVIEW ===\n")
cat("Sample size:", nrow(rdf), "\n")
cat("Number of areas:", length(unique(rdf$area)), "\n")  
cat("Number of schools:", length(unique(rdf$school)), "\n")
cat("Average students per school:", round(nrow(rdf)/length(unique(rdf$school)), 2), "\n")
cat("Average schools per area:", round(length(unique(rdf$school))/length(unique(rdf$area)), 2), "\n\n")

# =============================================================================
# FITTING MULTI-LEVEL NEGATIVE BINOMIAL MODELS (corresponding to big_model1, big_model2, big_model3)
# =============================================================================

# Multi-Level Model 1 (equivalent to big_model1 - no LBC factor)
cat("=== FITTING MULTI-LEVEL MODEL 1 (no LBC factor) ===\n")

big_mlm1 <- glmmTMB(delinquency ~ mgb + stegb + ltegb + 
                      vbs + nrwa + nrws + nrwp + pf + np + 
                      nle + sne + apc + ats + atp + D45 + D47 + 
                      bcv + ica + delinf + fses + age + 
                      gender + (1|area/school), 
                    data = rdf, 
                    family = nbinom2())

## Model summary
summary(big_mlm1)

## Get the rounded coefficients
round(fixef(big_mlm1)$cond, 3)

## Get the IRR and its confidence intervals
est1 <- get_irr_ci(big_mlm1)
print(round(est1, 3))

## Get ICC information
icc1 <- get_icc(big_mlm1)
cat("\nICC Area:", round(icc1$icc_area, 4))
cat("\nICC School:", round(icc1$icc_school, 4))
cat("\nICC Total:", round(icc1$icc_total, 4), "\n\n")

# Multi-Level Model 2 (equivalent to big_model2 - with binary LBC factor)
cat("=== FITTING MULTI-LEVEL MODEL 2 (with binary LBC factor) ===\n")

big_mlm2 <- glmmTMB(delinquency ~ mgb + stegb + ltegb + 
                      vbs + nrwa + nrws + nrwp + pf + np + 
                      nle + sne + apc + ats + atp + D45 + D47 + 
                      bcv + ica + delinf + fses + age + 
                      gender + lbc + (1|area/school), 
                    data = rdf, 
                    family = nbinom2())

## Model summary
summary(big_mlm2)

## Get the rounded coefficients
round(fixef(big_mlm2)$cond, 3)

## Get the IRR and its confidence intervals
est2 <- get_irr_ci(big_mlm2)
print(round(est2, 3))

## Get ICC information
icc2 <- get_icc(big_mlm2)
cat("\nICC Area:", round(icc2$icc_area, 4))
cat("\nICC School:", round(icc2$icc_school, 4))
cat("\nICC Total:", round(icc2$icc_total, 4), "\n\n")

## Likelihood Ratio Test (comparing mlm1 vs mlm2)
anova(big_mlm1, big_mlm2)

## Wald test for all fixed effects except intercept
# Extract fixed effects variance-covariance matrix
vcov_mlm2 <- vcov(big_mlm2)$cond
coef_mlm2 <- fixef(big_mlm2)$cond

# Test all coefficients except intercept (positions 2 to end)
n_coef <- length(coef_mlm2)
wald_stat2 <- as.numeric(t(coef_mlm2[2:n_coef]) %*% solve(vcov_mlm2[2:n_coef, 2:n_coef]) %*% coef_mlm2[2:n_coef])
wald_df2 <- n_coef - 1
wald_p2 <- 1 - pchisq(wald_stat2, wald_df2)

cat("Wald test for Model 2:\n")
cat("Chi-square =", round(wald_stat2, 3), ", df =", wald_df2, ", p =", round(wald_p2, 6), "\n\n")

# Multi-Level Model 3 (equivalent to big_model3 - with three-level LBC1 factor)
cat("=== FITTING MULTI-LEVEL MODEL 3 (with three-level LBC1 factor) ===\n")

big_mlm3 <- glmmTMB(delinquency ~ mgb + stegb + ltegb + 
                      vbs + nrwa + nrws + nrwp + pf + np + 
                      nle + sne + apc + ats + atp + D45 + D47 + 
                      bcv + ica + delinf + fses + age +
                      gender + lbc1 + (1|area/school), 
                    data = rdf, 
                    family = nbinom2())

## Model summary
summary(big_mlm3)

## Get the rounded coefficients
round(fixef(big_mlm3)$cond, 3)

## Get the IRR and its confidence intervals
est3 <- get_irr_ci(big_mlm3)
print(round(est3, 3))

## Get ICC information
icc3 <- get_icc(big_mlm3)
cat("\nICC Area:", round(icc3$icc_area, 4))
cat("\nICC School:", round(icc3$icc_school, 4))
cat("\nICC Total:", round(icc3$icc_total, 4), "\n\n")

## Likelihood Ratio Test (comparing mlm1 vs mlm3)
anova(big_mlm1, big_mlm3)

## Wald test for all fixed effects except intercept
vcov_mlm3 <- vcov(big_mlm3)$cond
coef_mlm3 <- fixef(big_mlm3)$cond

# Test all coefficients except intercept (positions 2 to end)
n_coef3 <- length(coef_mlm3)
wald_stat3 <- as.numeric(t(coef_mlm3[2:n_coef3]) %*% solve(vcov_mlm3[2:n_coef3, 2:n_coef3]) %*% coef_mlm3[2:n_coef3])
wald_df3 <- n_coef3 - 1
wald_p3 <- 1 - pchisq(wald_stat3, wald_df3)

cat("Wald test for Model 3:\n")
cat("Chi-square =", round(wald_stat3, 3), ", df =", wald_df3, ", p =", round(wald_p3, 6), "\n\n")

# Multi-Level Model 4 (built upon big_mlm3 - with LBC1 interactions)
cat("=== FITTING MULTI-LEVEL MODEL 4 (with LBC1 interactions) ===\n")

big_mlm4 <- glmmTMB(delinquency ~ mgb + stegb + ltegb + 
                      vbs + nrwa + nrws + nrwp + pf + np + 
                      nle + sne + apc + ats + atp + D45 + D47 + 
                      bcv + ica + delinf + fses + age +
                      gender + lbc1 + 
                      lbc1:mgb + lbc1:stegb + lbc1:ltegb + lbc1:vbs + 
                      lbc1:nrwa + lbc1:nrws + lbc1:nrwp + lbc1:pf + 
                      lbc1:np + lbc1:nle + lbc1:sne + (1|area/school), 
                    data = rdf, 
                    family = nbinom2(),
                    control = glmmTMBControl(optimizer = "nlminb", 
                                           optCtrl = list(iter.max = 1000, 
                                                         eval.max = 1000)))

## Model summary
summary(big_mlm4)

## Get the rounded coefficients
round(fixef(big_mlm4)$cond, 3)

## Get the IRR and its confidence intervals
est4 <- get_irr_ci(big_mlm4)
print(round(est4, 3))

## Get ICC information
icc4 <- get_icc(big_mlm4)
cat("\nICC Area:", round(icc4$icc_area, 4))
cat("\nICC School:", round(icc4$icc_school, 4))
cat("\nICC Total:", round(icc4$icc_total, 4), "\n\n")

## Likelihood Ratio Test (comparing mlm3 vs mlm4)
anova(big_mlm3, big_mlm4)

## Wald test for all fixed effects except intercept
vcov_mlm4 <- vcov(big_mlm4)$cond
coef_mlm4 <- fixef(big_mlm4)$cond

# Test all coefficients except intercept (positions 2 to end)
n_coef4 <- length(coef_mlm4)
wald_stat4 <- as.numeric(t(coef_mlm4[2:n_coef4]) %*% solve(vcov_mlm4[2:n_coef4, 2:n_coef4]) %*% coef_mlm4[2:n_coef4])
wald_df4 <- n_coef4 - 1
wald_p4 <- 1 - pchisq(wald_stat4, wald_df4)

cat("Wald test for Model 4:\n")
cat("Chi-square =", round(wald_stat4, 3), ", df =", wald_df4, ", p =", round(wald_p4, 6), "\n\n")

## Wald test specifically for interaction terms
# Extract positions of interaction terms (these come after main effects)
interaction_names <- names(coef_mlm4)[grep("lbc1:", names(coef_mlm4))]
interaction_positions <- which(names(coef_mlm4) %in% interaction_names)

if(length(interaction_positions) > 0) {
  wald_int_stat <- as.numeric(t(coef_mlm4[interaction_positions]) %*% 
                             solve(vcov_mlm4[interaction_positions, interaction_positions]) %*% 
                             coef_mlm4[interaction_positions])
  wald_int_df <- length(interaction_positions)
  wald_int_p <- 1 - pchisq(wald_int_stat, wald_int_df)
  
  cat("Wald test for LBC1 interaction terms:\n")
  cat("Chi-square =", round(wald_int_stat, 3), ", df =", wald_int_df, ", p =", round(wald_int_p, 6), "\n\n")
}

