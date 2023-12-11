# All the mechanic to obtain Zf using Sf_closed ----

#' Scaling term N(t, sigma) (also possible to express in N(lambda, sigma))
#'
#' @param t Current time (t > 0)
#' @param sigma positive parameter conveying variance information
#' @param type string such as "linear", "exponential"...
#' @returns Scaling term to apply in this case to obtain Zf from Sf
#' @export
scaling_term_func = function(t, sigma, type) {
  if(type == "linear") {
    return(t*sigma)
  }
}

#' Computation of S_lambda f_sigma(x) given (t, z), using Sf_closed_func
#'
#' @param t Current time (t > 0)
#' @param z Current normalized position in space (z in R)
#' @param sigma positive parameter conveying variance information
#' @param type string such as "linear", "exponential"...
#' @returns Retrieve the original known form using Sf_closed
#' @export
normalization_in_space_and_time_from_Sf_closed_func = function(t, z, sigma, type) {
  lambda = sigma/t # t = sigma/lambda; sigma = lambda*t
  x = lambda*z # within [-lambda/2, lambda/2) when z is in [-1/2, 1/2)
  Sf_closed_func(type, sigma, lambda)(x)
}

#' Obtain Zf_func from Sf_closed_func, i.e. a (t, z) --> Zf(t, z) mapping,
#' for each type and sigma
#'
#' @param type string such as "linear", "exponential"...
#' @param sigma positive parameter conveying variance information
#' @returns Obtain Zf_func using Sf_closed
#' @export
Zf_func_from_Sf_closed_func = function(type, sigma) {
  # Computation of Zf given (t, z), using Sf_closed_func
  # t Current time (t > 0)
  # z Current normalized position in space (z in R)
  function(t, z) {
    # normalization in time
    Zf1 = sapply(t, normalization_in_space_and_time_from_Sf_closed_func, z, sigma, type)
    # normalization in space
    lambda = sigma/t
    normalization = scaling_term_func(t, sigma, type)
    Zf2 = (Zf1 - 1/lambda)*normalization
    Zf2 = as.numeric(Zf2)
    return(Zf2)
  }
}

# The helpers used to retrieve the simplified forms of the expressions ----

#' For the Rectangular and Linear types, Delta condition as a function of t and z
#'
#' @param t Current time (t > 0)
#' @param z Current position (z in R)
#' @returns Whether this current value (t,z) fulfill the condition Delta,
#' used in the Rectangular and Linear cases
#' @export
Delta_func = function(t, z) {
  # abs(t-1/2)+abs(z) <= 1/2 # correct form, but sensitive to rounding
  round(abs(t-1/2)+abs(z)-1/2, 10) <= 0
}

# All the simplified forms of the expressions, without reference to Sf_closed ----

#' Obtain Zf_func directly, i.e. a (t, z) --> Zf(t, z) mapping,
#' for each type and sigma
#'
#' @param type string such as "linear", "exponential"...
#' @param sigma positive parameter conveying variance information
#' @returns Obtain Zf_func directly
#' @export
Zf_func = function(type, sigma) {
  # Direct computation of Zf given (t, z)
  # t Current time (t > 0)
  # z Current normalized position in space (z in R)
  if(type == "linear") {
    return(
      function(t, z) {
        # t-floor(t+1/2) is the shifted fractional part of t
        # t-floor(t) is the fractional part of t
        # z-floor(z+1/2) is the shifted fractional part of z
        -pmin((t-floor(t+1/2))^2, (t-floor(t))*(t-floor(t)-1) + abs(z-floor(z+1/2)))
      }
    )
  }
}
