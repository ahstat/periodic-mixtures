# rm(list = ls())
# devtools::document()
# devtools::test()
# devtools::build()
# library(periodicmixtures)

# library(testthat)
# library(dplyr)

N = 1000
lambda0 = runif(N, 0, 3)
sigma0 = runif(N, 0, 3)
x0 = runif(N, -3, 3)

for(k in 1:N) {
  print(k)
  lambda = lambda0[k]
  sigma = sigma0[k]
  x = x0[k]
  MM = floor(lambda/(2*pi*sigma))
  f1 = (1/lambda)+(2/lambda)*(cos((MM+1)*pi*x/lambda))*(sin(MM*pi*x/lambda))/(sin(pi*x/lambda))
  f2 = (1/lambda)*sin((2*MM+1)*pi*x/lambda)/sin(pi*x/lambda)
  testthat::expect_equal(f1, f2)
}

sigma = 1
lambda = seq(from = 0, to = 100, length.out = 1000)
MM = floor(lambda/(2*pi*sigma))
plot(lambda, (2*MM+1)*pi/lambda, type = "l", log = "y")





###



type = "sinc"
eps = 1e-5
N = 7
length.out = 1e5
K = 100L
verbose = TRUE
approx_is_Fourier = TRUE

sigmas_and_lambdas = sigmas_and_lambdas_func(N)
sigmas = sigmas_and_lambdas$sigmas
lambdas = sigmas_and_lambdas$lambdas
vals = rep(NA, N)
for(k in 1:N) {
  if(k %% 1 == 0) {
    if(verbose) {
      print(k)
    }
  }
  sigma = sigmas[k]
  lambda = lambdas[k]
  interval = seq(from = -K, to = K, by = 1L)

  elements = build_f_g_Ff_Fg(type, sigma)

  if(approx_is_Fourier) {
    Ff = elements$Ff
    Sf_approx = Sf_approx_Fourier_func(Ff, lambda, interval)
  } else {
    f = elements$f
    Sf_approx = Sf_approx_direct_func(f, lambda, interval)
  }

  Sf_closed = Sf_closed_func(type, sigma, lambda)

  # Interval depending of lambda the periodicity
  x = seq(from = -2*lambda, to = 1*lambda, length.out = length.out)
  if(verbose) {
    plot(x, Sf_approx(x), type = "l", main = k, ylim = range(c(Sf_closed(x), Sf_approx(x))))
    lines(x, Sf_closed(x), col = "red")
  }
  val = mean(abs(Sf_approx(x) - Sf_closed(x)))
  vals[k] = val
}
if(verbose) {
  print(max(vals))
}
expect_true(all(vals <= eps)) # difference between direct sum and closed form sum is small

