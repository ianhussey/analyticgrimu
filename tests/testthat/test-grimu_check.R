test_that("output has the documented two-part structure", {
  res <- grimu_check(n1 = 8, n2 = 8, p_reported = 0.017, digits = 3)
  expect_type(res, "list")
  expect_true(all(c("summary", "details") %in% names(res)))
  expect_s3_class(res$summary, "tbl_df")
  expect_equal(nrow(res$summary), 1)
})

test_that("an attainable p-value is flagged consistent", {
  res <- grimu_check(n1 = 8, n2 = 8, p_reported = 0.017, digits = 3)
  expect_true(res$summary$consistent)
  expect_true(res$summary$p_granularity_consistent)
})

test_that("a reported U inconsistent with the p-value is caught", {
  # p = 0.041 cannot be produced together with U = 30 at n1 = 10, n2 = 12
  res <- grimu_check(
    n1 = 10,
    n2 = 12,
    u_reported = 30,
    p_reported = 0.041,
    digits = 3
  )
  expect_false(res$summary$consistent)
  expect_false(res$summary$u_matches_p)
})

test_that("missing sample sizes raise an error", {
  expect_error(grimu_check(n1 = NA, n2 = 8, p_reported = 0.05, digits = 3))
})

test_that("a p-value outside [0, 1] raises an error", {
  expect_error(grimu_check(n1 = 8, n2 = 8, p_reported = 1.2, digits = 3))
})

test_that("digits must be supplied with a p-value", {
  expect_error(grimu_check(n1 = 8, n2 = 8, p_reported = 0.05, digits = NA))
})
