# synpdx

synpdx is an R package for nonlinear mixed-effects (NLME) modeling and interaction analysis in in vivo patient-derived xenograft (PDX) combination studies.  
It fits tumor growth and drug-effect models, computes subject-level interaction distances (Dᵢ), and performs chi-square and bootstrap hypothesis testing to classify synergism, additivity, or antagonism.

---

## Installation

You can install the development version of synpdx from GitHub using remotes:

```r
remotes::install_github("mothquitoman-kor/synpdx")

---

## Example

The package includes a curated example dataset:

```r
library(synpdx)

# Load example dataset
csv_path <- system.file(
  "extdata", "PDX_CRC_BYL719+binimetinib_curated.csv",
  package = "synpdx"
)
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Example run: auto model selection + random effects selection
res <- synpdx_fit(
  dat = df,
  control = "control",
  drug_a = "BYL719",
  drug_b = "binimetinib",
  combo = "combination",
  B = 100,
  seed = 2025,
  cutoff = "auto",
  model = "auto",
  select_random_effects = TRUE,
  out_dir = "synpdx_out"
)

# Inspect results
res$fit$selected_model_name   # chosen model
res$chisq$overall             # chi-square test result
res$bootstrap$summary         # bootstrap summary
