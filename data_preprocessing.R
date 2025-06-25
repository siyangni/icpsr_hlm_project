# ===============================================================================
# Data Preprocessing Script for ICP-H Project
# This script consolidates all variable creation and recoding steps
# from the multiple scripts in old_scripts/
# ===============================================================================

# Load required packages
library(dplyr)
library(readr)
library(tibble)
library(eeptools)
library(zoo)
library(Hmisc)

# Read original data
df <- read_csv("data/original.csv")

# ===============================================================================
# STEP 1: Create basic delinquency components (from df_rvadded.R)
# ===============================================================================

# Creating 15 new variables consisting the delinquency variable
df$b24j <- ifelse(df$B24j == 1, 0, 1)
df$b24k <- ifelse(df$B24k == 1, 0, 1)
df$b24l <- ifelse(df$B24l == 1, 0, 1)
df$b24m <- ifelse(df$B24m == 1, 0, 1)
df$b24n <- ifelse(df$B24n == 1, 0, 1)
df$b27a <- ifelse(df$B27a == 0, 0, 1)
df$b27b <- ifelse(df$B27b == 0, 0, 1)
df$b27c <- ifelse(df$B27c == 0, 0, 1)
df$b27d <- ifelse(df$B27d == 0, 0, 1)
df$b27e <- ifelse(df$B27e == 0, 0, 1)
df$b27f <- ifelse(df$B27f == 0, 0, 1)
df$b27g <- ifelse(df$B27g == 0, 0, 1)
df$b27h <- ifelse(df$B27h == 0, 0, 1)
df$b27i <- ifelse(df$B27i == 0, 0, 1)
df$b27j <- ifelse(df$B27j == 0, 0, 1)

# ===============================================================================
# STEP 2: Create goal-related variables (from df_add_failure_positive_goal.R)
# ===============================================================================

# Creating variable mgb (monetary goal blockage)
df <- df %>% 
  mutate(b28i = case_when(B28i == 0 ~ 0,
                          B28i == 1 | B28i == 2 ~ 1, 
                          B28i == 3 ~ 2, 
                          B28i == 4 ~ 3)) %>%
  mutate(b28j = case_when(B28j == 0 ~ 0,
                          B28j == 1 | B28j == 2 ~ 1, 
                          B28j == 3 ~ 2, 
                          B28j == 4 ~ 3)) %>%
  mutate(b28k = case_when(B28k == 0 ~ 0,
                          B28k == 1 | B28k == 2 ~ 1, 
                          B28k == 3 ~ 2, 
                          B28k == 4 ~ 3)) %>%
  mutate(b28r = case_when(B28r == 0 ~ 0,
                          B28r == 1 | B28r == 2 ~ 1, 
                          B28r == 3 ~ 2, 
                          B28r == 4 ~ 3)) %>%
  # Create mgb
  mutate(mgb = b28i + b28j + b28k + b28r)

# Create stegb (short term educational goal blockage)
df <- df %>% 
  mutate(stegb = case_when(B28d == 0 ~ 0,
                           B28d == 1 | B28d == 2 ~ 1, 
                           B28d == 3 ~ 2, 
                           B28d == 4 ~ 3))

# Create ltegb (long term educational goal blockage)
df <- df %>% 
  mutate(ltegb = D47 - D48) %>%
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
                           ltegb == 5 ~ 10))

# ===============================================================================
# STEP 3: Add negative confrontation variables (from add_con_negative.R)
# ===============================================================================

df <- df %>% 
  # victimization by strangers
  mutate(vbs = B25e + B25o + B25p + B25q + B25r + B25s) %>%
  # negative relations with adults
  mutate(nrwa = B25g + B25h + B25i + B25j + B25k + B25l + B25m + B25n) %>%
  # negative relations with siblings
  mutate(nrws = B25a + B25b + B25c + B25d) %>%  
  mutate(nrwp = case_when(B28q == 0 ~ 0,
                          B28q == 1 | B28q == 2 ~ 1,
                          B28q == 3 ~ 2,  
                          B28q == 4 ~ 3)) %>%  # negative relations with peers
  mutate(pf = B26a + B26b) # Parents fighting

# Collapsing neighborhood problem variable series
df <- df %>% add_column(f75a=NA, f75b=NA, f75c=NA, f75d=NA, 
                        f75e=NA, f75f=NA, f75g=NA)

for (i in list(grep("F75", colnames(df)))){
  for (x in list(grep("f75", colnames(df))))
    df[, x] <- ifelse(df[, i]==0, 0, 1)
}

# neighborhood problems
df <- df %>%
  mutate(np = f75a + f75b + f75c + f75d + f75e + f75f + f75g)

# ===============================================================================
# STEP 4: Add loss of positive stimuli variables (from add_loss_positive.R)
# ===============================================================================

# Function to create b-variables based on B-variables  
c_nle <- function(x){
  x[x==0] <- 0
  x[x==1|x==2] <- 1
  x[x==3] <- 2
  x[x==4] <- 3
  return(x)
}

# Create four variables constituting nle(negative life events)
df <- df %>%
  add_column(b28a=NA, b28b=NA, b28c=NA, b28h=NA, b28l=NA)

df$b28a <- c_nle(df$B28a)
df$b28b <- c_nle(df$B28b)
df$b28c <- c_nle(df$B28c)
df$b28h <- c_nle(df$B28h)
df$b28l <- c_nle(df$B28l)

# Creating nle
df <- df %>%
  mutate(nle = b28a + b28b + b28c + b28h + b28l)

# ===============================================================================
# STEP 5: Add negative emotion variables (from add_negative_emotion.R)
# ===============================================================================

# Function to create b-variables based on B-variables  
c_sne <- function(x){
  y <- vector()
  y[x==4] <- 0
  y[x==3] <- 1
  y[x==2] <- 2
  y[x==1] <- 3
  return(y)
}

# Create 6 new variables based on B-variables
df$b29a <- c_sne(df$B29a)
df$b29b <- c_sne(df$B29b)
df$b29c <- c_sne(df$B29c)
df$b29d <- c_sne(df$B29d)
df$b29e <- c_sne(df$B29e)
df$b29f <- c_sne(df$B29f)

# Create sne (situational negative emotions)
df <- df %>% 
  mutate(sne = b29a + b29b + b29c + b29d + b29e + b29f)

# ===============================================================================
# STEP 6: Add family socioeconomic status (from add_fses.R)
# ===============================================================================

# Adding fses (family socio-economic status)
df <- df %>%
  mutate(fses = A7a + A7b + A7c + A7d + A7e + A7f +
                A7g + A7h + A7i + A7j + A7k + A7l +
                A7m + A7n + A7o + A7p + A7q)

# ===============================================================================
# STEP 7: Add attachment variables (from add_attach_*.R scripts)
# ===============================================================================

# Function for attachment variable creation
c_apc <- function(x){
  y <- vector()
  y[x==1] <- 0
  y[x==2] <- 1
  y[x==3] <- 2
  y[x==4] <- 3
  return(y)
}

# Adding attachment to primary caregivers variables
df <- df %>%
  mutate(c33 = c_apc(C33), c34 = c_apc(C34), c35 = c_apc(C35), 
         c36 = c_apc(C36), c37 = c_apc(C37), c43 = c_apc(C43)) %>%
  # Create variable apc (attachment to primary caregivers)
  mutate(apc = c33 + c34 + c35 + c36 + c37 + c43)

# Adding attachment to school variables
df <- df %>%
  mutate(d49 = c_apc(D49), d50 = c_apc(D50), d51 = c_apc(D51), 
         d53 = c_apc(D53), d54 = c_apc(D54)) %>%
  # Create variable ats (attachment to school)
  mutate(ats = d49 + d50 + d51 + d53 + d54)

# Adding attachment to peers variables
df <- df %>%
  mutate(e62 = c_apc(E62), e63 = c_apc(E63), e65 = c_apc(E65)) %>%
  # Create variable atp (attachment to peers)
  mutate(atp = e62 + e63 + e65)

# ===============================================================================
# STEP 8: Add educational goal commitment (from ed_goal.R)
# ===============================================================================

# Function to create sceg (short term commitment to educational goal) based on D45  
d <- function(x){
  y <- vector()
  y[x==1] <- 4
  y[x==2] <- 3
  y[x==3] <- 2
  y[x==4] <- 1
  return(y)
}

# Adding sceg
df <- df %>% 
  mutate(sceg = d(D45))

# ===============================================================================
# STEP 9: Add conventional values and activities (from add_value_activity.R)
# ===============================================================================

# Function for creating ica (involvement in conventional activities)
ica <- function(x){
  y <- vector()
  y[x==1] <- 0
  y[x==2] <- 1
  y[x==3] <- 2
  y[x==4] <- 3
  return(y)
}

df <- df %>% 
  # Creating bcv (belief in conventional values)
  mutate(bcv = G81a + G81b + G81c) %>% 
  # Creating the six new vars that construct ica
  mutate(b22a = ica(B22a), b22b = ica(B22b), b22c = ica(B22c), 
         b22d = ica(B22d), b22e = ica(B22e), b22f = ica(B22f)) %>% 
  # Creating ica 
  mutate(ica = b22a + b22b + b22c + b22d + b22e + b22f)

# ===============================================================================
# STEP 10: Add delinquent friends variable (from add_delinf.R)
# ===============================================================================

# Function for creating variable delinf (delinquent friends)
d_friends <- function(x){
  y <- vector()
  y[x==3|x==0] <- 0
  y[x==1] <- 1
  y[x==2] <- 2
  return(y)
}

df <- df %>% 
  # creating the nine variables constructing delinf
  mutate(e72a=d_friends(E72a), e72b=d_friends(E72b), e72c=d_friends(E72c), 
         e72d=d_friends(E72d), e72e=d_friends(E72e), e72f=d_friends(E72f), 
         e72g=d_friends(E72g), e72h=d_friends(E72h), e72i=d_friends(E72i)) %>% 
  # Creating delinf
  mutate(delinf = e72a + e72b + e72c + e72d + e72e + 
                  e72f + e72g + e72h + e72i)

# ===============================================================================
# STEP 11: Add final delinquency variable (from bug_fix_add_delinquency.R)
# ===============================================================================

df <- df %>%
  mutate(delinquency = b24j + b24k + b24l + b24m + b24n + 
                       b27a + b27b + b27c + b27d + b27e + 
                       b27f + b27g + b27h + b27i + b27j)

# ===============================================================================
# STEP 12: Add left-behind children variable (from add_lbc.R)
# ===============================================================================

# Add a variable lbc where 0 stands for nlbc, 1 stands for lbc
df <- df %>% 
  mutate(lbc = ifelse(LBC1==1, 0, 1))

# Add a variable lbc1 where recode LBC1 from 1, 2, 3, 4 to 0, 1, 2, 3
df <- df %>% 
  mutate(lbc1 = case_when(LBC1 == 1 ~ 0,
                          LBC1 == 2 ~ 1,
                          LBC1 == 3 ~ 2,
                          LBC1 == 4 ~ 3))

# ===============================================================================
# STEP 13: Add age and gender variables (from add_age_gender.R)
# ===============================================================================

# Create variable gender
df <- df %>% 
  mutate(gender = ifelse(A1==1, 0, 1))

# Creating variable age based on A3
# Converting A3 to character vectors and assign it to a3
a3 <- as.character(df$A3)

# Extract birth year from A3
# For YYYYMM format, extract first 4 characters
# For YYYY format (2002, 2003), use the whole value
birth_year <- ifelse(nchar(a3) == 6, 
                     as.numeric(substr(a3, 1, 4)),
                     as.numeric(a3))

# Calculate age: questionnaire year (2017) minus birth year
age <- 2017 - birth_year

# Adding variable age to dataframe
df$age <- age

# ===============================================================================
# VALIDATION CHECKS
# ===============================================================================

# Print summary statistics for key variables
cat("=== VALIDATION CHECKS ===\n")
cat("Sample size:", nrow(df), "\n\n")

# Check key variables
variables_to_check <- c("gender", "age", "lbc", "delinquency", "apc", "ats", "atp", 
                       "fses", "sne", "nle", "delinf", "bcv", "ica", "sceg", 
                       "mgb", "stegb", "ltegb", "vbs", "nrwa", "nrws", "nrwp", "pf", "np")

for(var in variables_to_check) {
  if(var %in% names(df)) {
    cat("Variable:", var, "\n")
    print(table(df[[var]], useNA = "ifany"))
    cat("\n")
  }
}

# ===============================================================================
# Save the processed data
# ===============================================================================

write_csv(df, "data/processed.csv")