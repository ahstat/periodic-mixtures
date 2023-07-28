test_that("sinc outputs correct normalized value", {
  expect_equal(sinc(0), 1)
  # zeros are integers for the normalized function
  expect_equal(sinc(1:10), rep(0, 10))
  expect_equal(sinc(-10:-1), rep(0, 10))
})

test_that("coth outputs correct value", {
  x = seq(from = -5, to = 5, length.out = 380)
  expect_equal(coth(x), (exp(x)+exp(-x))/(exp(x)-exp(-x)))
})

test_that("phi outputs correct value for some known values", {
  # Values obtained from https://en.wikipedia.org/wiki/Theta_function
  # Quoting wikipedia:
  # "Proper credit for most of these results goes to Ramanujan. See Ramanujan's
  # lost notebook and a relevant reference at Euler function. The Ramanujan
  # results quoted at Euler function plus a few elementary operations give the
  # results below, so they are either in Ramanujan's lost notebook or follow
  # immediately from it. See also Yi (2004)".
  expect_equal(phi(exp(-pi)), pi^(1/4) / gamma(3/4))
  C = (pi^(1/4) / gamma(3/4))
  phi_1 = C
  phi_4 = C * (8^(1/4) + 2)/4
  phi_16 = C * (4 + 128^(1/4) + (1024*8^(1/4) + 1024*2^(1/4))^(1/4))/16
  expect_equal(phi(exp(-pi)), phi_1)
  expect_equal(phi(exp(-4*pi)), phi_4)
  expect_equal(phi(exp(-16*pi)), phi_16)
})

test_that("sigmas_and_lambdas_func outputs a list of numeric values", {
  N = 20
  out = sigmas_and_lambdas_func(N = 20, seed = 1234)
  expect_equal(names(out), c("sigmas", "lambdas"))
  expect_equal(length(out$sigmas), N)
  expect_equal(length(out$lambdas), N)
  expect_true(all(out$lambdas >= 0 & out$lambdas <= 10))
  expect_true(all(out$sigmas >= 0 & out$sigmas <= 10))
})
