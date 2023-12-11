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

test_that("tilde_func has some closed form expressions", {
  step = 1/2^0 # tested with 1/2^3

  # For Rectangular, Linear, Exponential types (x tilde)
  f = function(x, lambda, sigma) {
    tilde_func(x, lambda)
  }
  g = function(x, lambda, sigma) {
    x - lambda*floor(x/lambda + 1/2)
  }
  check_equal_func(f, g, step)

  # For Rectangular type (sigma breve)
  f = function(x, lambda, sigma) {
    sigma_tilde2 = tilde_func(sigma, 2*lambda)
    return(sigma_tilde2)
  }
  g = function(x, lambda, sigma) {
    sigma - 2*lambda*floor(sigma/(2*lambda) + 1/2)
  }
  check_equal_func(f, g, step)

  # For Linear type (sigma tilde)
  f = function(x, lambda, sigma) {
    sigma_tilde = tilde_func(sigma, lambda)
    return(sigma_tilde)
  }
  g = function(x, lambda, sigma) {
    sigma - lambda*floor(sigma/lambda + 1/2)
  }
  check_equal_func(f, g, step)
})
