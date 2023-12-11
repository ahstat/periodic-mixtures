test_that("sigmas_and_lambdas_func outputs a list of numeric values", {
  N = 20
  out = sigmas_and_lambdas_func(N = 20, seed = 1234)
  expect_equal(names(out), c("sigmas", "lambdas"))
  expect_equal(length(out$sigmas), N)
  expect_equal(length(out$lambdas), N)
  expect_true(all(out$lambdas >= 0 & out$lambdas <= 10))
  expect_true(all(out$sigmas >= 0 & out$sigmas <= 10))
})
