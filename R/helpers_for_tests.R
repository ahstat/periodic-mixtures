#' Sample vectors of sigma and lambda parameters, only used for tests
#'
#' @param N Length of (sigma,lambda) couples to generate
#' @param seed Random seed
#' @returns a list of `N` sigmas and of `N` lambdas took uniformly in [0,10]
sigmas_and_lambdas_func = function(N = 400, seed = 1234) {
  set.seed(seed)
  sigmas = runif(N, 0, 10)
  lambdas = runif(N, 0, 10)
  return(list(sigmas = sigmas, lambdas = lambdas))
}

#' Test of a generic formula, only used for tests of Sf
#'
#' @param formula A function to test
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f`
#' @param original_function ground truth function parametrized with `type`
#' @param boundary_step if NULL, the N_tested parameter is used. Otherwise,
#' create a grid with steps of `boundary_step`, to take e.g. 1, 1/2, 1/4, 1/8
#' to finer and finer granularity
#' @param N_tested if boundary_step is NULL, create random values for x, sigma,
#' sigma, and test the formula for `N_tested` repetitions
#' @param Sf_approx which approximation is chosen, use Sf_approx_direct_func for
#' linear and exponential since it converges quickly
#' @returns pass the test is the formula in valid in those cases
test_formula = function(formula, type = "linear", original_function = f_func,
                        boundary_step = NULL, N_tested = 1000,
                        Sf_approx = Sf_approx_direct_func) {
  if(is.null(boundary_step)) {
    N = N_tested
    set.seed(1)
    df = data.frame(x = runif(N, -3, 3),
                    lambda = abs(runif(N, -3, 3)),
                    sigma = abs(runif(N, -3, 3)))
  } else {
    step = boundary_step
    x = seq(from = -2, to = 2, by = step)
    lambda = seq(from = step, to = 2, by = step)
    sigma = seq(from = step, to = 2, by = step)
    df = expand.grid(x = x, lambda = lambda, sigma = sigma)
  }
  N = nrow(df)
  Sf_x_approx = rep(NA, N)
  Sf_x_formula = rep(NA, N)
  for(i in 1:N) {
    x = df$x[i]
    lambda = df$lambda[i]
    sigma = df$sigma[i]
    interval = seq(from=-10000L, to = 10000L, by = 1L) # can still be different from the closed formula for very small values
    Sf_x_approx[i] = Sf_approx(original_function(type, sigma), lambda, interval)(x)
    Sf_x_formula[i] = formula(x, sigma, lambda)

    if(N > 100) {
      print(i)
    }
    if(abs(Sf_x_approx[i] - Sf_x_formula[i]) > 1e-5) {
      print(paste(lambda, sigma, x))
    }
    expect_equal(Sf_x_approx[i], Sf_x_formula[i])
  }
}

#' Other test of a generic formula, only used for tests
#'
#' @param f predicted function, that should be a function of (x, lambda, sigma)
#' @param g ground truth function, that should be a function of (x, lambda, sigma)
#' @param step grid over x, lambda, sigma with steps `step`, to take e.g. 1, 1/2, 1/4, 1/8,
#' to finer and finer granularity
#' @param x_range range of the interval for x
#' @param lambda_max max of the interval for lambda
#' @param sigma_max max of the interval for sigma
#' @returns check the validity of f == g on the grid of values
check_equal_func = function(f, g, step = 1/2^3, x_range = c(-2, 2), lambda_max = 2, sigma_max = 2) {
  x = seq(from = x_range[1], to = x_range[2], by = step)
  lambda = seq(from = step, to = lambda_max, by = step)
  sigma = seq(from = step, to = sigma_max, by = step)
  df = expand.grid(x = x, lambda = lambda, sigma = sigma)
  N = nrow(df)
  for(i in 1:N) {
    if(N > 1000) {
      if(i %% 100 == 0) {
        print(paste0(i, "/", N))
      }
    }
    predicted = f(df$x[i], df$lambda[i], df$sigma[i])
    groundtruth = g(df$x[i], df$lambda[i], df$sigma[i])
    expect_equal(predicted, groundtruth)
  }
}

#' Other test of a generic formula with (t,z), only used for tests
#'
#' @param f predicted function, that should be a function of (t, z)
#' @param g ground truth function, that should be a function of (t, z)
#' @param step grid over t, z, and sigma with steps `step`, to take e.g. 1, 1/2, 1/4, 1/8,
#' to finer and finer granularity
#' @param z_range range of the interval for z
#' @param t_max max of the interval for lambda
#' @returns check the validity of f == g on the grid of values
#' @export
check_equal_tz_func = function(f, g, step = 1/2^3, z_range = c(-2, 2), t_max = 2) {
  t = seq(from = step, to = t_max, by = step)
  z = seq(from = z_range[1], to = z_range[2], by = step)
  df = expand.grid(z = z, t = t)
  N = nrow(df)
  for(i in 1:N) {
    if(N > 1000) {
      if(i %% 100 == 0) {
        print(paste0(i, "/", N))
      }
    }
    predicted = f(df$t[i], df$z[i])
    groundtruth = g(df$t[i], df$z[i])
    expect_equal(predicted, groundtruth)
  }
}
