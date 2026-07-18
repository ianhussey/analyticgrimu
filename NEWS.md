# analyticgrimu 0.2.1

* First public, CRAN-oriented release.
* `grimu_check()` triangulates the reported sample sizes, U statistic, and
  p-value, checking the p-value against every reporting convention the
  package models: exact doubled-tail (R / SciPy), 'SPSS' Exact / 'StatXact'
  deviation-from-mean, mid-p, and tie-corrected asymptotic approximations
  over every achievable tie-correction K value.
* `check_u_consistency()` checks a reported U statistic for bounds and
  granularity (integer / half-integer) consistency given the group sizes.
* `grimu_map_pvalues()` enumerates the attainable U statistics and their
  implied p-values under all modelled conventions.
* `grimu_saturation()` reports the p-value granularity (saturation) of a
  design as an estimate of the false-positive rate of a GRIM-U check.
* `plot_grimu_probability()` visualises the attainable p-value space for one
  or more sample-size pairs.
* Added package documentation, a vignette, a test suite, and continuous
  integration.
