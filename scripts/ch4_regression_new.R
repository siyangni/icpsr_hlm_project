# Negative Binomial regression
library(pacman)
p_load(tidyverse, Hmisc, MASS, sandwich, lmtest, coxphw, stargazer)

# A function res for calculating robust standard error.
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

# A function to paste to clipboard
tabout <- function(output){
  print(output)
  capture.output(output, file = "clipboard", append = FALSE, 
                 type = "output", split = FALSE)
  lines <- readClipboard()
  for(i in 1 : 5) {lines <- gsub("  ", " ", lines, fixed=TRUE)}
  lines <- gsub(" ", "\t", lines, fixed=TRUE)
  writeClipboard(lines)
}

# Load the data
load("data/imp_aggregate.RData")
df <- dat$imputations$imp1 # Take the first set of the five

# Preparing the operation data frame for regression
rdf <- df %>%
  dplyr::select(delinquency, mgb, stegb, ltegb, vbs, nrwa,
         nrws, nrwp, pf, np, nle, sne, apc,
         ats, atp, D45, D47, bcv, ica, delinf, fses,
         age, gender, school, lbc, lbc1)


# Checking Assumptions
## General summary
summary(rdf)

## Visualization
ggplot(rdf, aes(delinquency, fill = lbc1)) + geom_histogram(binwidth = 1) + 
  facet_grid(lbc1 ~ ., margins = TRUE, scales = "free")

## Conditional mean and variance table
with(rdf, tapply(delinquency, lbc1, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

## From the above table, we can see:
# 1. The means of delinquency vary by lbc1, suggesting 
# 2. Within each group, sd is larger than mean (overdispersion)


# Fitting Negative Binomial Model (whole sample with binary lbc)

summary(big_model1 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                               vbs + nrwa +nrws + nrwp + pf + np + 
                               nle + sne + apc + ats + atp + D45 + D47 + 
                               bcv + ica + delinf + fses + age + 
                               factor(gender), data = rdf))

summary(big_model2 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                         vbs + nrwa +nrws + nrwp + pf + np + 
                         nle + sne + apc + ats + atp + D45 + D47 + 
                         bcv + ica + delinf + fses + age + 
                         factor(gender) + factor(lbc),data = rdf))

## Get the rounded coefficients
round(coef(big_model2),3)

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(big_model2), confint(big_model2)))
print(est)

## Obtaining Cluster Robust standard error
crse_big2 <- coeftest(big_model2, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_big2[,2], 3)

## Likelihood Ratio Test
lrtest(big_model1, big_model2)

## Wald test 
wald(coef(big_model2),vcov(big_model2), 2:24)


# Fitting Negative Binomial Model (whole sample with three-level lbc1)

summary(big_model3 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                    vbs + nrwa +nrws + nrwp + pf + np + 
                    nle + sne + apc + ats + atp + D45 + D47 + 
                    bcv + ica + delinf + fses + age +
                    factor(gender) + factor(lbc1), data = rdf))

## Get the rounded coefficients
round(coef(big_model3),3)

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(big_model3), confint(big_model3)))
round(est[,1],3)

## Obtaining Cluster Robust standard error
crse_big3 <- coeftest(big_model3, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_big3[,2], 3)

## Likelihood Ratio Test
lrtest(big_model1, big_model3)

## Wald test 
wald(coef(big_model3),vcov(big_model3), 2:24)




# Fitting Negative Binomial Distribution (lbc two parents, model 1)

rdf2 <- rdf %>% filter(lbc1==2) %>% 
  dplyr::select(delinquency, mgb, stegb, ltegb, vbs, nrwa,
         nrws, nrwp, pf, np, nle, sne, apc,
         ats, atp, D45, D47, bcv, ica, delinf, fses,
         age, gender, school)

summary(s2m1 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                     vbs + nrwa + nrws + nrwp + pf + np + 
                     nle + apc + ats + atp + D45 + D47 + 
                     bcv + ica + delinf + fses + age + factor(gender), 
                     data = rdf2))

## Get the rounded coefficients
round(coef(s2m1),3)

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s2m1), confint(s2m1)))
round(est[,1],3)

## Obtaining Cluster Robust standard error
crse_s2m1 <- coeftest(s2m1, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s2m1[,2], 3)

## Wald test 
wald(coef(s2m1),vcov(s2m1), 2:22)


# Fitting Negative Binomial Distribution (lbc two parents, model 2)
summary(s2m2 <- glm.nb(delinquency ~ sne + apc + ats + atp + D45 + D47 + 
                       bcv + ica + delinf + fses + age + factor(gender), 
                       data = rdf2))

## Get the rounded coefficients
round(coef(s2m2),3)

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s2m2), confint(s2m2)))
round(est[,1],3)

## Obtaining Cluster Robust standard error
crse_s2m2 <- coeftest(s2m2, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s2m2[,2], 3)

## Wald test 
###wald.test(b = coef(s2m2), Sigma = vcov(s2m2), Terms = 2:12)
wald(coef(s2m2),vcov(s2m2), 2:13)


# Fitting Negative Binomial Distribution (lbc two parents, model 3)
summary(s2m3 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                       vbs + nrwa +nrws + nrwp + pf + np + 
                       nle + sne + apc + ats + atp + D45 + 
                       D47 + bcv + ica + delinf + fses + age +
                       factor(gender), data = rdf2))

## round results to three decimals 
round(summary(s2m3)$coefficients, 3) 

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s2m3), confint(s2m3)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s2m3 <- coeftest(s2m3, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s2m3[,2], 3)

## Wald test 
###wald.test(b = coef(s2m3), Sigma = vcov(s2m3), Terms = 2:23)
wald(coef(s2m3),vcov(s2m3), 2:23)


# Fitting Negative Binomial Distribution (lbc one parent, model 1)

rdf1 <- rdf %>% filter(lbc1==1) %>% 
  dplyr::select(delinquency, mgb, stegb, ltegb, vbs, nrwa,
                nrws, nrwp, pf, np, nle, sne, apc,
                ats, atp, D45, D47, bcv, ica, delinf, fses,
                age, gender, school)

summary(s1m1 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                         vbs + nrwa + nrws + nrwp + pf + np + 
                         nle + apc + ats + atp + D45 + D47 + 
                         bcv + ica + delinf + fses + age + factor(gender), 
                         data = rdf1))

## round results to three decimals 
round(summary(s1m1)$coefficients, 3) 

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s1m1), confint(s1m1)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s1m1 <- coeftest(s1m1, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s1m1[,2], 3)

## Wald test 
###wald.test(b = coef(s2m3), Sigma = vcov(s2m3), Terms = 2:23)
wald(coef(s1m1),vcov(s1m1), 2:22)


# Fitting Negative Binomial Distribution (lbc one parents, model 2)
summary(s1m2 <- glm.nb(delinquency ~ sne + apc + ats + atp + D45 + D47 + 
                       bcv + ica + delinf + fses + age + factor(gender), 
                       data = rdf1))

## round results to three decimals 
round(summary(s1m2)$coefficients, 3) 

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s1m2), confint(s1m2)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s1m2 <- coeftest(s1m2, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s1m2[,2], 3)

## Wald test 
wald(coef(s1m2),vcov(s1m2), 2:13)


# Fitting Negative Binomial Distribution (lbc one parent, model 3)

summary(s1m3 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                         vbs + nrwa +nrws + nrwp + pf + np + 
                         nle + sne + apc + ats + atp + D45 + 
                         D47 + bcv + ica + delinf + fses + age +
                         factor(gender), data = rdf1))

## round results to three decimals 
round(summary(s1m3)$coefficients, 3) 

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s1m3), confint(s1m3)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s1m3 <- coeftest(s1m3, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s1m3[,2], 3)

## Wald test 
wald(coef(s1m3),vcov(s1m3), 2:23)


# Fitting Negative Binomial Distribution (nlbc, model 1)

rdf0 <- rdf %>% filter(lbc1==0) %>% 
  dplyr::select(delinquency, mgb, stegb, ltegb, vbs, nrwa,
         nrws, nrwp, pf, np, nle, sne, apc,
         ats, atp, D45, D47, bcv, ica, delinf, fses,
         age, gender, school)

summary(s0m1 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                         vbs + nrwa +nrws + nrwp + pf + np + 
                         nle + apc + ats + atp + D45 + D47 + 
                         bcv + ica + delinf + fses + age + factor(gender), 
                         data = rdf0))

## round results to three decimals 
round(summary(s0m1)$coefficients, 3) 

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s0m1), confint(s0m1)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s0m1 <- coeftest(s0m1, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s0m1[,2], 3)

## Wald test 
wald(coef(s0m1),vcov(s0m1), 2:22)


# Fitting Negative Binomial Distribution (nlbc, model 2)
summary(s0m2 <- glm.nb(delinquency ~ sne + apc + ats + atp + D45 + D47 + 
                       bcv + ica + delinf + fses + age + factor(gender), 
                       data = rdf0))

## round results to three decimals 
round(summary(s0m2)$coefficients, 3) 

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s0m2), confint(s0m2)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s0m2 <- coeftest(s0m2, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s0m2[,2], 3)

## Wald test 
wald(coef(s0m2),vcov(s0m2), 2:13)


# Fitting Negative Binomial Distribution (nlbc, model 3)
summary(s0m3 <- glm.nb(delinquency ~ mgb + stegb + ltegb + 
                         vbs + nrwa + nrws + nrwp + pf + np + 
                         nle + sne + apc + ats + atp + D45 + 
                         D47 + bcv + ica + delinf + fses + age + 
                         factor(gender), data = rdf0))

## round results to three decimals 
round(summary(s0m3)$coefficients, 3)

## get the IRR and its confidence intervals
est <- exp(cbind(IRR = coef(s0m3), confint(s0m3)))
round(est,3)

## Obtaining Cluster Robust standard error
crse_s0m3 <- coeftest(s0m3, vcov. = vcovCL, cluster= ~school, type="HC0")
round(crse_s0m3[,2], 3)

## Wald test 
wald(coef(s0m3),vcov(s0m3), 2:23)
