# synpdx

synpdx is an R package for nonlinear mixed-effects (NLME) modeling and interaction analysis in *in vivo* patient-derived xenograft (PDX) combination studies.  
It fits tumor growth and drug-effect models, computes subject-level interaction distances (Dᵢ), and performs chi-square and bootstrap hypothesis testing to classify synergism, additivity, or antagonism.

---

## Installation

You can install the development version of synpdx from GitHub using remotes:

```r
remotes::install_github("mothquitoman-kor/synpdx")
```

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
```

---

## Output Files

When `out_dir` is specified, the following files are written:

| Output | Description |
|--------|-------------|
| `individual_prediction_plot.png` | Fitted vs observed tumor growth curves |
| `combo_pred_vs_obs.csv` | Predicted vs observed values for the combination arm |
| `Di_by_subject.csv` | Subject-level interaction distances (Dᵢ) |
| `chi_square_result.csv`, `chi_square_subject_results.csv`, `chi_square_test_plot.png` | Chi-square test outputs |
| `bootstrap_summary_result.csv`, `bootstrap_subject_Si.csv`, `bootstrap_ci_overlay.png` | Bootstrap test outputs |

---

## Citation

If you use synpdx, please cite:

> Something

