#' Normalized sinus cardinal function
#' @param x Input number or vector
#' @returns the numeric value sinc(x)
#' @export
sinc = function(x) {
  # function used for the Fourier transform
  out = sin(pi*x)/(pi*x)
  # replace any undefined value to 1
  if(sum(is.nan(out)) > 0) {
    out[is.nan(out)] = 1
  }
  return(out)
}

#' Hyperbolic cotangent function
#' @param x Input number or vector
#' @returns the numeric value coth(x)
#' @export
coth = function(x) {
  # function used for explicit values for polynomial
  return(1/tanh(x))
}

#' Phi function
#'
#' Phi function as defined in
#' https://en.wikipedia.org/wiki/Jacobi_theta_functions_(notational_variations)
#' and also known as the third theta function in
#' https://en.wikipedia.org/wiki/Theta_function
#' @param x Input number or vector
#' @param maxiter Maximum number of iterations used. Note that the series
#' generally converge very quickly
#' @returns the numeric value phi(x)
#' @export
phi = function(x, maxiter = 3e4) {
  # function used to get explicit values for the Gaussian type
  # theta3(z, q)
  elliptic::theta3(0, q = x, maxiter = maxiter)
}

#' Tilde function, to push any number x or sigma into [-lambda/2, lambda/2)
#'
#' @param x Input number (either x or sigma)
#' @param lambda lambda positive parameter corresponding to the periodicity
#' @returns x_tilde the number in [-lambda/2, lambda/2) such that
#' `x = x_tilde+i*lambda` (with i integer)
tilde_func = function(x, lambda) {
  tilde_inner_func = function(x, lambda) {
    x_tilde = x %% lambda
    if(x_tilde >= lambda/2) {
      x_tilde = x_tilde - lambda
    }
    return(x_tilde)
  }
  sapply(x, tilde_inner_func, lambda)
}

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

#' Test of a generic formula, only used for tests
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
  Sf_x_closed = rep(NA, N)
  Sf_x_formula = rep(NA, N)
  for(i in 1:N) {
    x = df$x[i]
    lambda = df$lambda[i]
    sigma = df$sigma[i]
    interval = interval=seq(from=-10000L, to = 10000L, by = 1L) # can still fail for very small values
    Sf_x_closed[i] = Sf_approx(original_function(type, sigma), lambda, interval)(x)
    Sf_x_formula[i] = formula(x, sigma, lambda)

    if(N > 100) {
      print(i)
    }
    expect_equal(Sf_x_closed[i], Sf_x_formula[i])
  }
}
