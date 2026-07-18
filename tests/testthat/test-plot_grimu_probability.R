test_that("a ggplot object is returned", {
  p <- plot_grimu_probability(n1_vec = c(8, 30), n2_vec = c(9, 31))
  expect_s3_class(p, "ggplot")
})

test_that("mismatched vector lengths raise an error", {
  expect_error(plot_grimu_probability(n1_vec = c(8, 30), n2_vec = 9))
})
