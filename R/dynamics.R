# All the mechanic to obtain Zf using Sf_closed ----

#' Scaling term N(t, sigma) (also possible to express in N(lambda, sigma))
#'
#' @param t Current time (t > 0)
#' @param sigma positive parameter conveying variance information
#' @param type string such as "linear", "exponential"...
#' @returns Scaling term to apply in this case to obtain Zf from Sf
#' @export
scaling_term_func = function(t, sigma, type) {
  if(type == "rectangular") {
    return(sigma)
  } else if(type == "linear") {
    return(t*sigma)
  } else if(type == "exponential") {
    return(t*sigma)
  } else if(type == "polynomial") {
    # normalization (exp(2*t)/(2*t))*sigma
    # == normalization `dominant term` (gives ~1/(2t) in 0)
    return((exp(2*t)/(2*t))*sigma)
  } else if(type == "polynomial_one_shift") {
    # normalization ((exp(2*t)-1)/(2*t))*sigma
    # == normalization `~1 in 0`
    return((exp(2*t)-1)/(2*t)*sigma)
  } else if(type == "polynomial_two_shifts") {
    # normalization ((exp(2*t)-1-2*t)/(2*t))*sigma
    # == normalization `~t in 0`
    return((exp(2*t)-1-2*t)/(2*t)*sigma)
  } else if(type == "gaussian") {
    # normalization (exp(pi*t^2)/(2*t))*sigma
    # == normalization `dominant term` (gives ~1/(2t) in 0) (yes, same as poly)
    return((exp(pi*t^2)/(2*t))*sigma)
  } else if(type == "gaussian_one_shift") {
    # normalization (exp(pi*t^2-1)/(2*t))*sigma
    # == normalization `at least ~t in 0`
    # (not only bounded; hence no need to check another normalization)
    return((exp(pi*t^2)-1)/(2*t)*sigma)
  } else if(type == "sinc") {
    return(sigma)
  } else if(type == "sinc2") {
    return(sigma)
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
    type0 = type
    type0 = ifelse(grepl("polynomial", type0), "polynomial", type0)
    type0 = ifelse(grepl("gaussian", type0), "gaussian", type0)
    Zf1 = sapply(t, normalization_in_space_and_time_from_Sf_closed_func, z, sigma, type0)
    # normalization in space
    lambda = sigma/t
    normalization = scaling_term_func(t, sigma, type)
    Zf2 = (Zf1 - 1/lambda)*normalization
    Zf2 = as.numeric(Zf2)
    return(Zf2)
  }
}

# The helpers used to retrieve the simplified forms of the expressions ----

#' For the Linear type, Delta condition as a function of t and z
#'
#' @param t Current time (t > 0)
#' @param z Current position (z in R)
#' @returns Whether this current value (t,z) fulfill the condition Delta,
#' used in the Linear case
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
  if(type == "rectangular") {
    return(
      function(t, z) {
        # z-floor(z+1/2) is the shifted fractional part of z
        frac_t_over_2 = t/2-floor(t/2) # == {t/2}
        value_diff = round(abs(frac_t_over_2-1/2)+abs(z-floor(z+1/2))-1/2, 10)
        conditionB = (frac_t_over_2 < 1/2)
        conditionA = (value_diff<0) | ((value_diff==0) & (conditionB))
        (1-2*frac_t_over_2)  + (1-conditionA)*(-1)^(conditionB)
      }
    )
  } else if(type == "linear") {
    return(
      function(t, z) {
        # t-floor(t+1/2) is the shifted fractional part of t
        # t-floor(t) is the fractional part of t
        # z-floor(z+1/2) is the shifted fractional part of z
        -pmin((t-floor(t+1/2))^2, (t-floor(t))*(t-floor(t)-1) + abs(z-floor(z+1/2)))
      }
    )
  } else if(type == "exponential") {
    return(
      function(t, z) {
        res = t*cosh((1-2*abs(z-floor(z+1/2)))/t)/sinh(1/t)-t^2
        idx=which(t==0) # case multiple t and single z
        if(length(idx) > 0) {
          if(length(t) == 1) { # case with single t and multiple z
            idx=1:length(res)
          }
          res[idx] = 0
        }
        res
      }
    )
  } else if(type == "polynomial") {
    return(
      function(t, z) {
        # normalization (exp(2*t)/(2*t))*sigma
        # == normalization `dominant term` (gives ~1/(2t) in 0, not convenient)
        (cos(2*pi*z)*exp(2*t) - 1)/(exp(2*t)+exp(-2*t)-2*cos(2*pi*z))
      }
    )
  } else if(type == "polynomial_one_shift") {
    return(
      function(t, z) {
        # normalization ((exp(2*t)-1)/(2*t))*sigma
        # == normalization `~1 in 0`
        (1 - 1/exp(2*t))*(cos(2*pi*z)*exp(2*t) - 1)/(exp(2*t)+exp(-2*t)-2*cos(2*pi*z))
      }
    )
  } else if(type == "polynomial_two_shifts") {
    return(
      function(t, z) {
        # normalization ((exp(2*t)-1-2*t)/(2*t))*sigma
        # == normalization `~t in 0`
        (1-exp(-2*t)*(1+2*t))*(cos(2*pi*z)*exp(2*t) - 1)/(exp(2*t)-2*cos(2*pi*z)+exp(-2*t))
      }
    )
  } else if(type == "gaussian") {
    return(function(t, z) {
      # version 3:
      sapply(t, function(t){Re((exp(pi*t^2)/2)*(jacobi::jtheta3(z = pi*z, q = exp(-pi*t^2))-1))})
      #       # normalization (exp(pi*t^2)/(2*t))*sigma
      #       # Zf = cos(2*pi*z) + sum_{k=2}^{+inf} exp(-pi*t^2*(k^2-1))*cos(2*k*pi*z)
      #       maxiter = 3e4
      #
      #       # Debug code:
      #       # z=0
      #       # t=seq(from = 0.3, to = 5, length.out = 1000)[-1]
      #
      #       # Version 1: Not stable from t = 3:
      #       # out = theta3(z = pi*z, q = exp(-pi*t^2), maxiter = maxiter)
      #       # v1 = (out - 1)*exp(pi*t^2)/2
      #
      #       # Version 2: Stable even for t > 3
      #       out = cos(2*pi*z)
      #       out_previous = -2
      #
      #       q = exp(-pi*t^2)
      #       for (k in 2:maxiter) {
      #         out_new = out + q^(k^2-1) * cos(2*k*pi*z)
      #         if(elliptic::near.match(out, out_new) & k >= 5 & elliptic::near.match(out, out_previous)) {
      #           return(out)
      #         }
      #
      #         out_previous = out
      #         out = out_new
      #       }
      #       return(NA) # maximum iterations reached
      #       # Debug code: compare v1 and v2 over time
      #       # v2 = out
      #       # plot(t, v1, type = "l")
      #       # lines(t, v2, col = "red")
      #       # plot(t, v1-v2, type = "l")
    })
  } else if(type == "gaussian_one_shift") {
    return(function(t, z) {
      sapply(t, function(t){(1/2)*(exp(pi*t^2)-1)*Re(jacobi::jtheta3(z = pi*z, q = exp(-pi*t^2))-1)})
    })
  } else if(type == "sinc") {
    return(
      function(t, z) {
        # normalization sigma

        # first formula (complex in 0)
        #out = sapply(t, function(t){t*(sin(pi*z*(2*floor(1/(2*t))+1))/sin(pi*z) - 1)})
        #if(sum(is.nan(out)) > 0) {
        #  out[is.nan(out)] = (2*t*floor(1/(2*t)))[is.nan(out)]
        #}
        #return(out)

        # but there is easier, from the proof of Sf for sinc:
        # S - 1/lambda = (2/lambda)*sum_{k=1}^{floor(lambda/(2sigma))} cos(2 pi k x/lambda)
        #              = (1/sigma)*2*t*sum_{k=1}^{floor(1/(2t))} cos(2 k pi z)
        # so with N=sigma:
        # Zf = 2*t*sum_{k=1}^{floor(1/(2t))} cos(2 k pi z)
        # to go t=0.0001 we need 5000 terms
        sapply(t, function(t){if(t>1/2) {
            return(0)
          } else {
            return(sapply(z, function(z){2*t*sum(cos(2*(1:floor(1/(2*t)))*pi*z))})) # normalization sigma
            #return(sapply(z, function(z){2*sum(cos(2*(1:floor(1/(2*t)))*pi*z))})) # normalization lambda discarded, not good shape
          }})

        # The derivative w.r.t t is: dZf/dt = 0 with:
        #   sum_{k=1}^{k_max} cos(2 k pi z) = 0
      }
    )
  } else if(type == "sinc2") {
    function(t, z) {
    sapply(t, function(t){if(t > 1) {
      return(0)
    } else {
      N <- floor(1/t)
      return(sapply(z, function(z){2*t*sum((1 - (1:N)*t) * cos(2*(1:N)*pi*z))}))
    }})
    }
  } else {
    stop("Zf closed-form is not defined for this type in Zf_func")
  }
}
