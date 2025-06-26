# Commence multiple imputation via Amelia II 
library(pacman)
p_load(tidyverse, dplyr, Amelia, funModeling)

df <- read_csv("data/processed.csv")

df <- df %>%
  dplyr::select( b24j, b24k, b24l, b24m, b24n, b27a, b27b, 
          b27c, b27d, b27e, b27f, b27g, b27h, b27i, 
          b27j, b28i, b28j, b28k, b28r, 
          stegb, D47, D48, 
          B25e, B25o, B25p, B25q, B25r, B25s,  
          B25g, B25h, B25i, B25j, B25k, B25l, B25m, B25n,
          B25a, B25b, B25c, B25d, 
          nrwp, B26a, B26b,                    
          f75a, f75b, f75c, f75d, f75e, f75f, f75g, 
          b28a, b28b, b28c, b28h, b28l,
          b29a, b29b, b29c, b29d, b29e, b29f,
          A7a, A7b, A7c, A7d, A7e, A7f, A7g, A7h, A7i, 
          A7j, A7k, A7l, A7m, A7n, A7o, A7p, A7q,
          c33, c34, c35, c36, c37, c43,
          d49, d50, d51, d53, d54,
          e62, e63, e65, D45, D47, 
          G81a, G81b, G81c,
          b22a, b22b, b22c, b22d, b22e, b22f,
          e72a, e72b, e72c, e72d, e72e, e72f, e72g, e72h, e72i,
          LBC1, ID, Ir, area, school, age, gender)

# Check NA values again and see if all columns are numeric
status(df)
# Change df from tbl to dataframe so Amelia can process
df <- as.data.frame(df)

# Setting logical bound (value has to be positive)
bds <- matrix(data=0, nrow = 117, ncol = 3)

for (i in 1:117){
  bds[i,1]=i
  bds[i,2]=0
  bds[i,3]=Inf
}

# Generate 5 imputed datasets
complete_data <-amelia(df, m=5, cs='area', p2s=1, 
                       # treating vars constructing delinquency as ordinal
                       # treating D47, D48 as ordinal 
                       ords=c("b24j", "b24k", "b24l", "b24m",
                              "b24n", "b27a", "b27b", "b27c",
                              "b27d", "b27e", "b27f", "b27g", 
                              "b27h", "b27i", "b27j", "D47", "D48"),
                       # treating LBC1 and gender as nominal 
                       noms=c("LBC1", "gender"),
                       # excluding ID, investigator ID and school
                       idvars=c("ID", "Ir", "school"),
                       # Imposing logical bound 
                       bounds = bds) 

# Output five separate datasets
write.amelia(obj=complete_data, separate = TRUE, 
             file.stem= "data/imputed", impvar = "imp", 
             orig.data = TRUE)

# Output a single dataset
# write.amelia(obj=complete_data , separate = FALSE, 
            # file.stem= "data/imputed_all", impvar
            # = "imp", orig.data = TRUE)

# Save amelia output object
save(complete_data, file = "data/imputations.RData")




# Construct aggregate variables for the research
dat <- complete_data
## dataset 1
dat$imputations$imp1 <- dat$imputations$imp1 %>% # adding delinquency
  mutate(delinquency = b24j + b24k + b24l + b24m + b24n + 
           b27a + b27b + b27c + b27d + b27e + 
           b27f + b27g + b27h + b27i + b27j,
         # monetary goal blockage
         mgb = b28i + b28j + b28k + b28r,
         # long-term educational goal blockage
         ltegb = D47 - D48
  ) %>% 
  mutate(ltegb = case_when(ltegb == -5 ~ 0,
                           ltegb == -4 ~ 1,
                           ltegb == -3 ~ 2,  
                           ltegb == -2 ~ 3,
                           ltegb == -1 ~ 4,
                           ltegb == 0 ~ 5,
                           ltegb == 1 ~ 6,
                           ltegb == 2 ~ 7,
                           ltegb == 3 ~ 8,
                           ltegb == 4 ~ 9,
                           ltegb == 5 ~ 10)) %>% 
  # victimization by strangers
  mutate(vbs = B25e + B25o + B25p + B25q + B25r + B25s,
         # negative relations with adults
         nrwa = B25g + B25h + B25i + B25j + B25k + B25l + B25m + B25n,
         # negative relations with siblings
         nrws = B25a + B25b + B25c + B25d,
         # Parents fighting
         pf = B26a + B26b,
         # neighborhood problems
         np = f75a + f75b + f75c + f75d + f75e + f75f + f75g,
         # negative life events
         nle = b28a + b28b + b28c + b28h + b28l,
         # situational negative emotion
         sne = b29a + b29b + b29c + b29d + b29e + b29f,
         # family socio-economic status
         fses = A7a + A7b + A7c + A7d + A7e + A7f +
           A7g + A7h + A7i + A7j + A7k + A7l +
           A7m + A7n + A7o + A7p + A7q,
         # attachment to primary caregiver
         apc = c33 + c34 + c35 + c36 + c37 + c43,
         # attachment to school
         ats = d49 + d50 + d51 + d53 + d54,
         # attachment to peers
         atp = e62 + e63 + e65,
         # belief in conventional values
         bcv = G81a + G81b + G81c,
         # involvement in conventional activities
         ica = b22a + b22b + b22c + b22d + b22e + b22f,
         # delinquent friends 
         delinf = e72a + e72b + e72c + e72d + e72e + 
           e72f + e72g + e72h + e72i,
         # lbc (non lbc or lbc)
         lbc = ifelse(LBC1==1, 0, 1),
         # lbc1 (non lbc, lbc1l, lbc2)
         lbc1=case_when(LBC1==1 ~ 0,
                      LBC1==2 | LBC1==3 ~ 1,
                      LBC1==4 ~ 2))


## dataset 2
dat$imputations$imp2 <- dat$imputations$imp2 %>% # adding delinquency
  mutate(delinquency = b24j + b24k + b24l + b24m + b24n + 
           b27a + b27b + b27c + b27d + b27e + 
           b27f + b27g + b27h + b27i + b27j,
         # monetary goal blockage
         mgb = b28i + b28j + b28k + b28r,
         # long-term educational goal blockage
         ltegb = D47 - D48
  ) %>% 
  mutate(ltegb = case_when(ltegb == -5 ~ 0,
                           ltegb == -4 ~ 1,
                           ltegb == -3 ~ 2,  
                           ltegb == -2 ~ 3,
                           ltegb == -1 ~ 4,
                           ltegb == 0 ~ 5,
                           ltegb == 1 ~ 6,
                           ltegb == 2 ~ 7,
                           ltegb == 3 ~ 8,
                           ltegb == 4 ~ 9,
                           ltegb == 5 ~ 10)) %>% 
  # victimization by strangers
  mutate(vbs = B25e + B25o + B25p + B25q + B25r + B25s,
         # negative relations with adults
         nrwa = B25g + B25h + B25i + B25j + B25k + B25l + B25m + B25n,
         # negative relations with siblings
         nrws = B25a + B25b + B25c + B25d,
         # Parents fighting
         pf = B26a + B26b,
         # neighborhood problems
         np = f75a + f75b + f75c + f75d + f75e + f75f + f75g,
         # negative life events
         nle = b28a + b28b + b28c + b28h + b28l,
         # situational negative emotion
         sne = b29a + b29b + b29c + b29d + b29e + b29f,
         # family socio-economic status
         fses = A7a + A7b + A7c + A7d + A7e + A7f +
           A7g + A7h + A7i + A7j + A7k + A7l +
           A7m + A7n + A7o + A7p + A7q,
         # attachment to primary caregiver
         apc = c33 + c34 + c35 + c36 + c37 + c43,
         # attachment to school
         ats = d49 + d50 + d51 + d53 + d54,
         # attachment to peers
         atp = e62 + e63 + e65,
         # belief in conventional values
         bcv = G81a + G81b + G81c,
         # involvement in conventional activities
         ica = b22a + b22b + b22c + b22d + b22e + b22f,
         # delinquent friends 
         delinf = e72a + e72b + e72c + e72d + e72e + 
           e72f + e72g + e72h + e72i,
         # lbc (non lbc or lbc)
         lbc = ifelse(LBC1==1, 0, 1),
         # lbc1 (non lbc, lbc1l, lbc2)
         lbc1=case_when(LBC1==1 ~ 0,
                        LBC1==2 | LBC1==3 ~ 1,
                        LBC1==4 ~ 2))

## Dataset 3
dat$imputations$imp3 <- dat$imputations$imp3 %>% # adding delinquency
  mutate(delinquency = b24j + b24k + b24l + b24m + b24n + 
           b27a + b27b + b27c + b27d + b27e + 
           b27f + b27g + b27h + b27i + b27j,
         # monetary goal blockage
         mgb = b28i + b28j + b28k + b28r,
         # long-term educational goal blockage
         ltegb = D47 - D48
  ) %>% 
  mutate(ltegb = case_when(ltegb == -5 ~ 0,
                           ltegb == -4 ~ 1,
                           ltegb == -3 ~ 2,  
                           ltegb == -2 ~ 3,
                           ltegb == -1 ~ 4,
                           ltegb == 0 ~ 5,
                           ltegb == 1 ~ 6,
                           ltegb == 2 ~ 7,
                           ltegb == 3 ~ 8,
                           ltegb == 4 ~ 9,
                           ltegb == 5 ~ 10)) %>% 
  # victimization by strangers
  mutate(vbs = B25e + B25o + B25p + B25q + B25r + B25s,
         # negative relations with adults
         nrwa = B25g + B25h + B25i + B25j + B25k + B25l + B25m + B25n,
         # negative relations with siblings
         nrws = B25a + B25b + B25c + B25d,
         # Parents fighting
         pf = B26a + B26b,
         # neighborhood problems
         np = f75a + f75b + f75c + f75d + f75e + f75f + f75g,
         # negative life events
         nle = b28a + b28b + b28c + b28h + b28l,
         # situational negative emotion
         sne = b29a + b29b + b29c + b29d + b29e + b29f,
         # family socio-economic status
         fses = A7a + A7b + A7c + A7d + A7e + A7f +
           A7g + A7h + A7i + A7j + A7k + A7l +
           A7m + A7n + A7o + A7p + A7q,
         # attachment to primary caregiver
         apc = c33 + c34 + c35 + c36 + c37 + c43,
         # attachment to school
         ats = d49 + d50 + d51 + d53 + d54,
         # attachment to peers
         atp = e62 + e63 + e65,
         # belief in conventional values
         bcv = G81a + G81b + G81c,
         # involvement in conventional activities
         ica = b22a + b22b + b22c + b22d + b22e + b22f,
         # delinquent friends 
         delinf = e72a + e72b + e72c + e72d + e72e + 
           e72f + e72g + e72h + e72i,
         # lbc (non lbc or lbc)
         lbc = ifelse(LBC1==1, 0, 1),
         # lbc1 (non lbc, lbc1l, lbc2)
         lbc1=case_when(LBC1==1 ~ 0,
                        LBC1==2 | LBC1==3 ~ 1,
                        LBC1==4 ~ 2))

## dataset 4
dat$imputations$imp4 <- dat$imputations$imp4 %>% # adding delinquency
  mutate(delinquency = b24j + b24k + b24l + b24m + b24n + 
           b27a + b27b + b27c + b27d + b27e + 
           b27f + b27g + b27h + b27i + b27j,
         # monetary goal blockage
         mgb = b28i + b28j + b28k + b28r,
         # long-term educational goal blockage
         ltegb = D47 - D48
  ) %>% 
  mutate(ltegb = case_when(ltegb == -5 ~ 0,
                           ltegb == -4 ~ 1,
                           ltegb == -3 ~ 2,  
                           ltegb == -2 ~ 3,
                           ltegb == -1 ~ 4,
                           ltegb == 0 ~ 5,
                           ltegb == 1 ~ 6,
                           ltegb == 2 ~ 7,
                           ltegb == 3 ~ 8,
                           ltegb == 4 ~ 9,
                           ltegb == 5 ~ 10)) %>% 
  # victimization by strangers
  mutate(vbs = B25e + B25o + B25p + B25q + B25r + B25s,
         # negative relations with adults
         nrwa = B25g + B25h + B25i + B25j + B25k + B25l + B25m + B25n,
         # negative relations with siblings
         nrws = B25a + B25b + B25c + B25d,
         # Parents fighting
         pf = B26a + B26b,
         # neighborhood problems
         np = f75a + f75b + f75c + f75d + f75e + f75f + f75g,
         # negative life events
         nle = b28a + b28b + b28c + b28h + b28l,
         # situational negative emotion
         sne = b29a + b29b + b29c + b29d + b29e + b29f,
         # family socio-economic status
         fses = A7a + A7b + A7c + A7d + A7e + A7f +
           A7g + A7h + A7i + A7j + A7k + A7l +
           A7m + A7n + A7o + A7p + A7q,
         # attachment to primary caregiver
         apc = c33 + c34 + c35 + c36 + c37 + c43,
         # attachment to school
         ats = d49 + d50 + d51 + d53 + d54,
         # attachment to peers
         atp = e62 + e63 + e65,
         # belief in conventional values
         bcv = G81a + G81b + G81c,
         # involvement in conventional activities
         ica = b22a + b22b + b22c + b22d + b22e + b22f,
         # delinquent friends 
         delinf = e72a + e72b + e72c + e72d + e72e + 
           e72f + e72g + e72h + e72i,
         # lbc (non lbc or lbc)
         lbc = ifelse(LBC1==1, 0, 1),
         # lbc1 (non lbc, lbc1l, lbc2)
         lbc1=case_when(LBC1==1 ~ 0,
                        LBC1==2 | LBC1==3 ~ 1,
                        LBC1==4 ~ 2))

## Dataset 5
dat$imputations$imp5 <- dat$imputations$imp5 %>% # adding delinquency
  mutate(delinquency = b24j + b24k + b24l + b24m + b24n + 
           b27a + b27b + b27c + b27d + b27e + 
           b27f + b27g + b27h + b27i + b27j,
         # monetary goal blockage
         mgb = b28i + b28j + b28k + b28r,
         # long-term educational goal blockage
         ltegb = D47 - D48
  ) %>% 
  mutate(ltegb = case_when(ltegb == -5 ~ 0,
                           ltegb == -4 ~ 1,
                           ltegb == -3 ~ 2,  
                           ltegb == -2 ~ 3,
                           ltegb == -1 ~ 4,
                           ltegb == 0 ~ 5,
                           ltegb == 1 ~ 6,
                           ltegb == 2 ~ 7,
                           ltegb == 3 ~ 8,
                           ltegb == 4 ~ 9,
                           ltegb == 5 ~ 10)) %>% 
  # victimization by strangers
  mutate(vbs = B25e + B25o + B25p + B25q + B25r + B25s,
         # negative relations with adults
         nrwa = B25g + B25h + B25i + B25j + B25k + B25l + B25m + B25n,
         # negative relations with siblings
         nrws = B25a + B25b + B25c + B25d,
         # Parents fighting
         pf = B26a + B26b,
         # neighborhood problems
         np = f75a + f75b + f75c + f75d + f75e + f75f + f75g,
         # negative life events
         nle = b28a + b28b + b28c + b28h + b28l,
         # situational negative emotion
         sne = b29a + b29b + b29c + b29d + b29e + b29f,
         # family socio-economic status
         fses = A7a + A7b + A7c + A7d + A7e + A7f +
           A7g + A7h + A7i + A7j + A7k + A7l +
           A7m + A7n + A7o + A7p + A7q,
         # attachment to primary caregiver
         apc = c33 + c34 + c35 + c36 + c37 + c43,
         # attachment to school
         ats = d49 + d50 + d51 + d53 + d54,
         # attachment to peers
         atp = e62 + e63 + e65,
         # belief in conventional values
         bcv = G81a + G81b + G81c,
         # involvement in conventional activities
         ica = b22a + b22b + b22c + b22d + b22e + b22f,
         # delinquent friends 
         delinf = e72a + e72b + e72c + e72d + e72e + 
           e72f + e72g + e72h + e72i,
         # lbc (non lbc or lbc)
         lbc = ifelse(LBC1==1, 0, 1),
         # lbc1 (non lbc, lbc1l, lbc2)
         lbc1=case_when(LBC1==1 ~ 0,
                        LBC1==2 | LBC1==3 ~ 1,
                        LBC1==4 ~ 2))

# Output
save(dat, file = "data/imp_aggregate.RData")
