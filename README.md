# Bayesian Modelling of Combined Obesity Prevalence (Sex-Stratified Framework)

![R](https://img.shields.io/badge/R-4.0%2B-blue)
![Stan](https://img.shields.io/badge/Stan-Enabled-brightgreen)
![Reproducibility](https://img.shields.io/badge/Reproducible-Yes-success)

---

## Overview
This repository implements a **Bayesian hierarchical modelling framework** for estimating and integrating sex-stratified prevalence, with specific application to **combined obesity estimation** in epidemiological studies.

The workflow supports:
- Integration of male and female prevalence estimates
- Hierarchical modelling across sampling locations
- Robust uncertainty quantification using Bayesian inference

This pipeline underpins published work on obesity prevalence among Ghanaian populations.

---

## Associated Publication
This pipeline is associated with the following study:

**Obirikorang, C., Adu, E.A., Anto, E.O. et al. (2024)**  
*Prevalence and risk factors of obesity among undergraduate student population in Ghana: an evaluation study of body composition indices.*  
BMC Public Health 24, 877.  
https://doi.org/10.1186/s12889-023-17175-5

---

## Scientific Context
Accurate estimation of obesity prevalence requires:
- Integration of heterogeneous datasets (e.g., sex-stratified cohorts)
- Adjustment for site-level variability
- Explicit modelling of uncertainty

This framework addresses these challenges using a **Bayesian hierarchical model**, enabling:
- Partial pooling across study sites
- Improved estimates under sparse or imbalanced data
- Joint inference across demographic strata

---

## Repository Structure
```
.
├── combine_model_female.R      # Female-specific Bayesian model execution
├── combine_model_male.R        # Male-specific Bayesian model execution
├── combine_prevalence.stan     # Core hierarchical Bayesian model
├── data/                       # Input datasets (user-provided)
├── results/                    # Model outputs and summaries
└── README.md
```

---

## Data Requirements

### Input Data
Tabular prevalence or obesity classification data stratified by sex and site.

### Required Fields
| Column Name       | Description                          |
|------------------|--------------------------------------|
| location / site   | Sampling location identifier         |
| sample_size       | Number of individuals sampled        |
| positive_cases    | Number classified as obese           |
| sex               | Male or Female                       |

---

## Model Specification

The model is implemented in `combine_prevalence.stan`.

### Likelihood
y_i ~ Binomial(n_i, p_i)

### Linear Predictor
logit(p_i) = alpha + u_site[i]

### Random Effects
u_site ~ Normal(0, sigma)

Where:
- y_i = number of obese individuals  
- n_i = total sample size  
- p_i = prevalence probability  
- u_site = site-level variation  

---

## Installation

### Prerequisites
- R ≥ 4.0
- Stan (rstan or cmdstanr recommended)

### Install Dependencies
```r
install.packages(c("rstan", "tidyverse", "posterior"))
```

---

## Usage

### Run Female Model
```r
source("combine_model_female.R")
```

### Run Male Model
```r
source("combine_model_male.R")
```

---

## Outputs

The pipeline generates:
- Posterior prevalence estimates
- 95% credible intervals
- Convergence diagnostics (R-hat, ESS)

Stored in:
```
results/
```

---

## Reproducibility

- Fully script-driven workflow  
- Explicit probabilistic model (Stan)  
- Deterministic with fixed seeds  
- Suitable for replication and extension  

---

## Extensibility

The framework can be extended to:
- Multi-level models (region → site)
- Inclusion of covariates (age, lifestyle, diet)
- Longitudinal or temporal modelling

---

## Citation

If you use this pipeline, please cite:

Obirikorang, C., Adu, E.A., Anto, E.O. et al. (2024).  
Prevalence and risk factors of obesity among undergraduate student population in Ghana.  
BMC Public Health 24, 877.  
https://doi.org/10.1186/s12889-023-17175-5

---

## Contact
For questions or collaboration, please open an issue in the repository.

