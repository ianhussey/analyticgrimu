test_that("output has the expected structure", {
  res <- grimu_map_pvalues(n1 = 8, n2 = 9)
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("U", "is_integer", "p_exact") %in% names(res)))
  expect_true(any(grepl("^p_corr_K", names(res))))
  expect_true(any(grepl("^p_uncorr_K", names(res))))
})

test_that("U values lie on the half-integer lattice within bounds", {
  res <- grimu_map_pvalues(n1 = 8, n2 = 9, u_min = 0, u_max = 72)
  expect_true(all(res$U >= 0 & res$U <= 72))
  expect_true(all(abs(res$U * 2 - round(res$U * 2)) < 1e-9))
})

test_that("all p-value columns stay within [0, 1]", {
  res <- grimu_map_pvalues(n1 = 8, n2 = 9)
  p_cols <- grep("^p_", names(res), value = TRUE)
  for (col in p_cols) {
    vals <- res[[col]]
    vals <- vals[!is.na(vals)]
    expect_true(all(vals >= 0 & vals <= 1), info = paste("column", col))
  }
})

test_that("exact p-values are NA for half-integer U", {
  res <- grimu_map_pvalues(n1 = 8, n2 = 9)
  expect_true(all(is.na(res$p_exact[!res$is_integer])))
})

test_that("empty input sizes return an empty tibble", {
  res <- grimu_map_pvalues(n1 = NA, n2 = 9)
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 0)
})
