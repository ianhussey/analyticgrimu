# analyticgrimu: Granularity Testing for Mann-Whitney U Test Results

<!-- badges: start -->
[![R-CMD-check](https://github.com/ianhussey/analyticgrimu/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ianhussey/analyticgrimu/actions/workflows/R-CMD-check.yaml)

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![DOI](https://zenodo.org/badge/1131604883.svg)](https://doi.org/10.5281/zenodo.21434617)
<!-- badges: end -->

`analyticgrimu` is a trustworthiness-assessment / forensic-metascience tool for deciding whether a reported Mann-Whitney U test result is mathematically possible. Given the two group sizes, it recalculates the attainable U statistics and p-values and compares them against the reported values, taking into account (a) the rounding of the reported inputs and (b) the several p-value conventions used by different statistical software.

It implements the GRIM-U method of Heathers and Grimes ([Medical Evidence Project, Report 1, 2025](https://medicalevidenceproject.org/grim-u-observation-establish-impossible-p-values-ranked-tests/)), which extends the GRIM family of granularity tests from means to ranked tests. Grimes' original implementation is available [here](https://github.com/drg85/GRIMU).

## The problem

For two groups of sizes `n1` and `n2`, the Mann-Whitney U statistic can take only a discrete set of values: integers when there are no ties, and half-integers when there are. Each attainable U maps, in turn, to a discrete set of attainable p-values. Reported statistics are rounded, so each true value lies within a small window around the printed number.

Question: does *any* attainable value fall inside the rounding window of the reported one?

| answer | verdict |
|---|---|
| no attainable U / p / (U, p) pair fits the window | **`inconsistent`** (the reported result cannot have arisen from this design) |
| some attainable value fits | **`consistent`** |

An `inconsistent` verdict means the reported statistic, p-value, or their pairing is impossible for the stated sample sizes under every reporting convention the package models.

## The method

`grimu_check()` performs a triangulation across three quantities:

1. **N &harr; U** — is the reported U within `[0, n1*n2]` and on the integer / half-integer lattice?
2. **N &harr; p** — is the reported p-value attainable by *some* U for this design?
3. **U &harr; p** — can the *specific* reported U produce the *specific* reported p-value?

The candidate p-value set checked at each U covers the conventions used across common software:

- **exact, doubled smaller tail** (R `wilcox.test()`, SciPy);
- **exact, deviation-from-mean** (SPSS Exact, StatXact);
- **mid-p**;
- **asymptotic** (continuity-corrected and uncorrected) under every achievable tie-correction `K = sum(t_i^3 - t_i)`, matching the tie-corrected variance used by R, SPSS, Stata, and SciPy.

For one-sided alternatives both group labellings are enumerated and OR'd, because which group is labelled "1" determines whether the reported one-sided p is the left- or right-tail value and this convention differs between packages.

### What an inconsistency does and does not prove

An `inconsistent` verdict means the reported values cannot be produced by the modelled Mann-Whitney U procedures at the stated sample sizes and rounding. Benign explanations should be ruled out before any integrity inference:

- **exact permutation tests with ties.** Tools such as `coin::wilcox_test(distribution = "exact")`, SAS `PROC NPAR1WAY EXACT`, StatXact, and the SPSS Exact Tests add-on compute the exact permutation p-value *conditional on the observed tie pattern*, which is not recoverable from `(n1, n2)`. A flag against a report from one of these tools may reflect this coverage gap rather than a real error.
- **transcription and rounding of the reported values** beyond the modelled window.

Inconsistent results can also have more problematic causes, including calculation errors and research-integrity issues. Treat a flag as a prompt to check the source, not as a verdict on its own.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ianhussey/analyticgrimu")
```

## Usage

See also the vignette in the R package (`vignette("analyticgrimu")`).

```r
library(analyticgrimu)

# Check a reported p-value against a design.
res <- grimu_check(n1 = 8, n2 = 8, p_reported = 0.017, digits = 3)
res$summary[, c("n1", "n2", "p_reported", "consistent")]
#> # A tibble: 1 x 4
#>      n1    n2 p_reported consistent
#>   <dbl> <dbl>      <dbl> <lgl>
#> 1     8     8      0.017 TRUE

# Triangulate a reported U and p-value together.
res_u <- grimu_check(n1 = 10, n2 = 12, u_reported = 30,
                     p_reported = 0.041, digits = 3)
res_u$summary[, c("consistent", "u_matches_p")]
#> # A tibble: 1 x 2
#>   consistent u_matches_p
#>   <lgl>      <lgl>
#> 1 FALSE      FALSE
# U = 30 and p = 0.041 are each attainable, but not together.
```

### Checking a U statistic on its own

```r
check_u_consistency(n1 = 10, n2 = 12, u_reported = 30)
#> # A tibble: 1 x 3
#>   u_bounds_consistent u_granularity_consistent u_possible
#>   <lgl>               <lgl>                    <lgl>
#> 1 TRUE                TRUE                     TRUE
```

### The attainable p-value space

```r
# Enumerate attainable U values and their implied p-values.
grimu_map_pvalues(n1 = 8, n2 = 9)

# Granularity (saturation) of a design: the proportion of possible rounded
# p-values the test can actually attain. Low values -> low false-positive rate.
grimu_saturation(n1 = 8, n2 = 9, decimals = 3)

# Visualise the attainable p-value space for several designs.
plot_grimu_probability(n1_vec = c(8, 30), n2_vec = c(9, 31))
```

## API

| function | purpose |
|---|---|
| `grimu_check(n1, n2, u_reported = NA, p_reported, digits, ...)` | main entry point; triangulates N, U, and p and returns a `summary` / `details` list |
| `check_u_consistency(n1, n2, u_reported)` | bounds and granularity check for a reported U statistic |
| `grimu_map_pvalues(n1, n2, u_min = NULL, u_max = NULL, alternative = "two.sided")` | enumerate attainable U values and their p-values under every modelled convention |
| `grimu_saturation(n1, n2, decimals = 3, ...)` | proportion of possible rounded p-values a design can attain (an estimate of the GRIM-U false-positive rate) |
| `plot_grimu_probability(n1_vec, n2_vec)` | visualise the attainable p-value space for one or more designs |

## References

- J. Heathers and D. R. Grimes (2025), "The GRIM-U observation: establishing impossible p-values in ranked tests," *Medical Evidence Project*, Report 1. <https://medicalevidenceproject.org/grim-u-observation-establish-impossible-p-values-ranked-tests/>
- N. J. L. Brown and J. A. J. Heathers (2017), "The GRIM test: a simple technique detects numerous anomalies in the reporting of results in psychology," *Social Psychological and Personality Science*, 8(4):363–369. [doi:10.1177/1948550616673876](https://doi.org/10.1177/1948550616673876)

## Suggested citation

Hussey, I. (2026). analyticgrimu: Granularity Testing for Mann-Whitney U Test Results. [Computer software] <https://github.com/ianhussey/analyticgrimu> [doi:10.5281/zenodo.21434617](https://doi.org/10.5281/zenodo.21434617)

## License

Code is MIT licensed © Ian Hussey (2026).

Text and images are CC BY 4.0 (see the suggested citation above).
