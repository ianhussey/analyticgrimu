test_that("saturation is a single proportion in [0, 1]", {
  s <- grimu_saturation(n1 = 8, n2 = 9, decimals = 3)
  expect_length(s, 1)
  expect_type(s, "double")
  expect_gte(s, 0)
  expect_lte(s, 1)
})

test_that("coarser rounding gives higher saturation", {
  # Fewer bins to fill -> a larger share of them is attainable.
  s_coarse <- grimu_saturation(n1 = 8, n2 = 9, decimals = 2)
  s_fine <- grimu_saturation(n1 = 8, n2 = 9, decimals = 3)
  expect_gte(s_coarse, s_fine)
})
