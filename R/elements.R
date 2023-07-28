#' Build f, g, Ff, and Fg functions according to the `f_g_Ff_Fg.R` file
#' using pre-defined types
#'
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f` and all the other functions
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type. The value is defined such that the Fourier tranform
#' sums to one if this parameter is set to `NULL`
#' @returns a list with the functions f, g, Ff, Fg (function, derivative and
#' the corresponding Fourier transforms), along with the type name as a string
#' and the sigma > 0 value
#' @export
build_f_g_Ff_Fg = function(type, sigma=NULL) {
  if(is.null(sigma)) {
    sigma = sigma_such_as_Fourier_tranform_sums_to_one_func(type)
  }
  f = f_func(type, sigma)
  g = g_func(type, sigma)
  Ff = Ff_func(type, sigma)
  Fg = Fg_func(type, sigma)
  elements = list(type=type, sigma=sigma, f=f, g=g, Ff=Ff, Fg=Fg)
  return(elements)
}

#' Build f, g, Ff, and Fg functions according to the `f_g_Ff_Fg.R` file
#' from scratch
#'
#' @param f function f taking any real `x` as input and giving `f(x)`
#' @param g the derivative of f
#' @param Ff the Fourier transform of f
#' @param Fg the Fourier transform of g, can be set to NULL if Ff is defined
#' @param type optional string giving a name to the function f
#' @param sigma optional parameter for this specific f (f specified before
#' should not be sigma dependent)
#' @returns a list with the functions f, g, Ff, Fg (function, derivative and
#' the corresponding Fourier transforms), along optional type and sigma
#' @export
build_f_g_Ff_Fg_from_scratch = function(f=NULL, g=NULL, Ff=NULL, Fg=NULL, type=NULL, sigma=NULL) {
  elements = list(type=type, sigma=sigma, f=f, g=g, Ff=Ff, Fg=Fg)
  if(is.null(Fg) & !is.null(Ff)) {
    Fg = function(xi) {
      2 * pi * 1i * xi * Ff(xi)
    }
    elements$Fg = Fg
  }
  return(elements)
}

#' Build Sf, Sg, functions according to the `Sf_Sg.R` file, using built `elements`
#'
#' @param elements a list with the functions f, g, Ff, Fg (function, derivative and
#' the corresponding Fourier transforms), along with type and sigma. Some elements
#' of the list can be NULL
#' @param lambda positive parameter conveying the periodicity
#' @param position_interval terms k of the wrapped sum, ideally from -Inf to +Inf, but
#' approximated on the finite subset `position_interval`
#' @param frequency_interval terms k of the Poisson summation formula, ideally from -Inf to +Inf, but
#' approximated on the finite subset `frequency_interval`. Here we assume that `frequency_interval` is of form
#' [-K, K] for a certain integer `K`, otherwise it throws an error
#' @returns a list with the existing terms, possibly NULL (f, g, Ff, Fg, type, sigma)
#' with additional elements lambda, position_interval, frequency_interval,
#' and functions (possibly NULL):
#'     Sf_approx_direct, Sf_approx_Fourier, Sf_closed,
#' and Sg_approx_direct, Sg_approx_Fourier, Sg_closed.
#' @export
build_Sf_Sg = function(elements,
                       lambda = 1,
                       position_interval = seq(from = -5L, to = 5L, by = 1L),
                       frequency_interval = position_interval) {
  f = elements$f
  g = elements$g
  Ff = elements$Ff
  type = elements$type
  sigma = elements$sigma

  elements$lambda = lambda
  elements$position_interval = position_interval
  elements$frequency_interval = frequency_interval

  # Sf and Sg (approx direct)
  if(!is.null(f) & !is.null(position_interval)) {
    elements$Sf_approx_direct = Sf_approx_direct_func(f, lambda, position_interval)
  } else {
    elements$Sf_approx_direct = NULL
  }
  if(!is.null(g) & !is.null(position_interval)) {
    elements$Sg_approx_direct = Sg_approx_direct_func(g, lambda, position_interval)
  } else {
    elements$Sg_approx_direct = NULL
  }

  # Sf and Sg (approx Fourier)
  if(!is.null(Ff) & !is.null(frequency_interval)) {
    elements$Sf_approx_Fourier = Sf_approx_Fourier_func(Ff, lambda, frequency_interval)
    elements$Sg_approx_Fourier = Sg_approx_Fourier_func(Ff, lambda, frequency_interval)
  } else {
    elements$Sf_approx_Fourier = NULL
    elements$Sg_approx_Fourier = NULL
  }

  # Sf and Sg (closed form if known/available)
  if(!is.null(type)) {
    elements$Sf_closed = Sf_closed_func(type, sigma, lambda)
    elements$Sg_closed = Sg_closed_func(type, sigma, lambda)
  } else {
    elements$Sf_closed = NULL
    elements$Sg_closed = NULL
  }
  return(elements)
}
