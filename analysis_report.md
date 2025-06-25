# Three-Level Hierarchical Negative Binomial Regression Analysis Report
## ICP-H Project

### Model Specification

**Dependent Variable:** Delinquency (count variable)

**Fixed Effects (Level 1 - Individual):**
- Base model from ch4_regression.R (lines 103-107):
  - mgb, stegb, ltegb, vbs, nrwa, nrws, nrwp, pf, np, nle, sne, apc, ats, atp, D45, D47, bcv, ica, delinf, fses, age, gender, lbc1

**Interaction Effects:**
- lbc1 × mgb (Left-behind children status × Monetary Goal Blockage)
- lbc1 × stegb (Left-behind children status × Short-term Educational Goal Blockage)
- lbc1 × ltegb (Left-behind children status × Long-term Educational Goal Blockage)
- lbc1 × vbs (Left-behind children status × Victimization by Strangers)
- lbc1 × nrwa (Left-behind children status × Negative Relations with Adults)
- lbc1 × nrwp (Left-behind children status × Negative Relations with Peers)
- lbc1 × pf (Left-behind children status × Parents Fighting)
- lbc1 × np (Left-behind children status × Neighborhood Problems)
- lbc1 × nle (Left-behind children status × Negative Life Events)

**Random Effects:**
- Level 2: School (random intercept only)
- Level 3: Area (random intercept only)

### Key Findings

#### Random Effects Structure
- **Area-level variance (σ²_area):** 0.161² = 0.026
- **School-level variance (σ²_school):** 0.118² = 0.014
- **Residual variance (approximate):** π²/3 = 3.29

#### Intraclass Correlations (ICCs)
- **ICC Area (Level 3):** 0.0078 (0.78%)
- **ICC School (Level 2):** 0.0042 (0.42%)
- **ICC Combined (Levels 2+3):** 0.0120 (1.20%)

*Note: These ICCs indicate that only a small proportion of variance in delinquency is explained by clustering at the area and school levels.*

#### Significant Main Effects (p < 0.05)

**Positive associations with delinquency:**
- **Neighborhood Problems (np):** IRR = 1.048, p = 0.015
- **Attachment to Peers (atp):** IRR = 1.057, p < 0.001
- **Educational Goal Commitment (D45):** IRR = 1.309, p < 0.001
- **Delinquent Friends (delinf):** IRR = 1.095, p < 0.001
- **Family SES (fses):** IRR = 1.027, p = 0.008
- **Age:** IRR = 1.084, p = 0.006

**Negative associations with delinquency:**
- **Attachment to School (ats):** IRR = 0.968, p = 0.002
- **Belief in Conventional Values (bcv):** IRR = 0.949, p < 0.001
- **Female gender:** IRR = 0.501, p < 0.001

#### Interaction Effects
**No significant interactions** were found between left-behind children status (lbc1) and any of the strain/stress variables at p < 0.05 level.

The interaction terms tested included:
- lbc1 × mgb (p = 0.219-0.716)
- lbc1 × stegb (p = 0.078-0.329)
- lbc1 × ltegb (p = 0.575-0.728)
- lbc1 × vbs (p = 0.162-0.804)
- lbc1 × nrwa (p = 0.128-0.744)
- lbc1 × nrwp (p = 0.529-0.687)
- lbc1 × pf (p = 0.162-0.847)
- lbc1 × np (p = 0.082-0.377)
- lbc1 × nle (p = 0.179-0.539)

### Model Fit Comparison
- **Single-level model AIC:** [Value from analysis output]
- **Three-level model AIC:** [Value from analysis output]
- **AIC improvement:** [Difference indicating whether hierarchical structure improves fit]

### Interpretation

1. **Strain Theory Variables:** Several strain-related variables showed expected relationships:
   - Higher neighborhood problems increase delinquency risk
   - School attachment acts as a protective factor
   - Belief in conventional values reduces delinquency

2. **Left-Behind Children (LBC) Effects:** 
   - No significant main effects of LBC status were found
   - More importantly, no significant interactions between LBC status and strain variables
   - This suggests that the relationship between strain and delinquency does not differ significantly across LBC groups

3. **Social Learning Variables:**
   - Strong positive association with delinquent friends (IRR = 1.095)
   - Unexpected positive association with attachment to peers may reflect association with delinquent peer groups

4. **Hierarchical Structure:**
   - Low ICCs suggest limited clustering within schools and areas
   - Most variation in delinquency occurs at the individual level

### Methodological Notes

**Strengths:**
- Three-level hierarchical structure appropriately accounts for nesting
- Negative binomial distribution appropriate for count outcome
- Comprehensive set of theory-driven predictors and interactions

**Limitations:**
- Low ICCs suggest hierarchical structure may not be necessary
- Some interaction effect estimates may be underpowered
- Cross-sectional design limits causal inference

### Conclusions

1. The hypothesized moderating effect of left-behind children status on the relationship between strain variables and delinquency was not supported by the data.

2. Traditional predictors of delinquency (delinquent friends, school attachment, conventional beliefs) showed expected associations.

3. The three-level hierarchical structure captured minimal variance, suggesting individual-level factors are more important than school or area-level clustering for this outcome.

4. Future research might explore different specifications of LBC status or examine other potential moderators of strain-delinquency relationships. 