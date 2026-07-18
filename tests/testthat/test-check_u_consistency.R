test_that("a valid integer U passes both checks", {
  res <- check_u_consistency(n1 = 10, n2 = 12, u_reported = 30)
  expect_s3_class(res, "tbl_df")
  expect_true(res$u_bounds_consistent)
  expect_true(res$u_granularity_consistent)
  expect_true(res$u_possible)
})

test_that("a half-integer U is granularity-consistent (ties allowed)", {
  res <- check_u_consistency(n1 = 10, n2 = 12, u_reported = 30.5)
  expect_true(res$u_granularity_consistent)
  expect_true(res$u_possible)
})

test_that("a non-half-integer U fails the granularity check", {
  res <- check_u_consistency(n1 = 10, n2 = 12, u_reported = 30.3)
  expect_false(res$u_granularity_consistent)
  expect_false(res$u_possible)
})

test_that("U above n1 * n2 is out of bounds", {
  res <- check_u_consistency(n1 = 10, n2 = 12, u_reported = 200)
  expect_false(res$u_bounds_consistent)
  expect_false(res$u_possible)
})

test_that("a negative U is out of bounds", {
  res <- check_u_consistency(n1 = 10, n2 = 12, u_reported = -1)
  expect_false(res$u_bounds_consistent)
  expect_false(res$u_possible)
})
