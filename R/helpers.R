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
