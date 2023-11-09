#' Sum of f (direct sum approximation)
#'
#' @param f function f
#' @param lambda positive parameter conveying the periodicity
#' @param interval terms k of the wrapped sum, ideally from -Inf to +Inf, but
#' approximated on the finite subset `interval`
#' @returns sum of f taking a real `x` input and giving `Sf(x)` for this `f` and `lambda`
#' approximated from the direct sum on `interval`
#' @export
Sf_approx_direct_func = function(f, lambda, interval = seq(from = -5L, to = 5L, by = 1L)) {
  sum_function_single_x = function(x) {
    each_term = f(x + interval * lambda) # each term of the interval
    sum(each_term)
  }
  function(x) {
    sapply(x, sum_function_single_x)
  }
}

#' Sum of g (direct sum approximation)
#'
#' @param g function g
#' @param lambda positive parameter conveying the periodicity
#' @param interval terms k of the wrapped sum, ideally from -Inf to +Inf, but
#' approximated on the finite subset `interval`
#' @returns sum of g taking a real `x` input and giving `Sg(x)` for this `g` and `lambda`
#' approximated from the direct sum on `interval`
#' @export
Sg_approx_direct_func = Sf_approx_direct_func

#' Sum of f (Fourier sum approximation)
#'
#' @param Ff Fourier transform of the function f
#' @param lambda positive parameter conveying the periodicity
#' @param interval terms k of the Poisson summation formula, ideally from -Inf to +Inf, but
#' approximated on the finite subset `interval`. Here we assume that `interval` is of form
#' [-K, K] for a certain integer `K`, otherwise it throws an error
#' @returns sum of f taking a real `x` input and giving `Sf(x)` for this `f` and `lambda`
#' approximated from the Fourier sum on `interval`
#' @export
Sf_approx_Fourier_func = function(Ff, lambda, interval = seq(from = -5L, to = 5L, by = 1L)) {
  if(!0 %in% interval) {
    stop("The interval should be of form (-K):K for Sf_approx_Fourier_func")
  }
  if(length(interval) %% 2 != 1) {
    stop("The interval should be of form (-K):K for Sf_approx_Fourier_func")
  }
  K = (length(interval)-1)/2
  if(!all.equal(interval, (-K):K)) {
    stop("The interval should be of form (-K):K for Sf_approx_Fourier_func")
  }
  if(K > 0) {
    interval_positive = 1:K
  } else {
    interval_positive = c()
  }

  sum_function_single_x = function(x) {
    first_term = (1/lambda) * Ff(0)
    common_term = (2/lambda) * Ff(interval_positive/lambda)
    cos_in_x_term = cos(2*pi*interval_positive*x/lambda)
    each_term = c(first_term, common_term * cos_in_x_term)
    sum(each_term)
  }
  function(x) {
    sapply(x, sum_function_single_x)
  }
}

#' Sum of g (Fourier sum approximation)
#'
#' @param Ff Fourier transform of the function f (Ff is used even for computing here Sg)
#' @param lambda positive parameter conveying the periodicity
#' @param interval terms k of the Poisson summation formula, ideally from -Inf to +Inf, but
#' approximated on the finite subset `interval`. Here we assume that `interval` is of form
#' [-K, K] for a certain integer `K`, otherwise it throws an error
#' @returns sum of g taking a real `x` input and giving `Sg(x)` for this `g` and `lambda`
#' approximated from the Fourier sum on `interval`
#' @export
Sg_approx_Fourier_func = function(Ff, lambda, interval = seq(from = -5L, to = 5L, by = 1L)) {
  if(!0 %in% interval) {
    stop("The interval should be of form (-K):K for Sg_approx_Fourier_func")
  }
  if(length(interval) %% 2 != 1) {
    stop("The interval should be of form (-K):K for Sg_approx_Fourier_func")
  }
  K = (length(interval)-1)/2
  if(!all.equal(interval, (-K):K)) {
    stop("The interval should be of form (-K):K for Sg_approx_Fourier_func")
  }
  if(K > 0) {
    interval_positive = 1:K
  } else {
    interval_positive = c()
  }

  sum_function_single_x = function(x) {
    first_term = 0 # assuming the density function is odd and real
    common_term = -((4*pi)/lambda^2) * Ff(interval_positive/lambda)
    sin_in_x_term = interval_positive * sin(2*pi*interval_positive*x/lambda)
    each_term = c(first_term, common_term * sin_in_x_term)
    sum(each_term)
  }
  function(x) {
    sapply(x, sum_function_single_x)
  }
}

#' Sum of f (closed form)
#'
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc"
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @param lambda positive parameter conveying the periodicity
#' @returns sum of f taking a real `x` input and giving `Sf(x)` for this
#' `(type, sigma, lambda)` on its closed-form, or `NULL` if the closed-form
#' does not exist or is unknown
#' @export
Sf_closed_func = function(type, sigma, lambda) {
  if(is.null(type)) {
    return(NULL)
  } else if(type == "rectangular") {
    sum_function = function(x) {
      x = tilde_func(x, lambda)
      sigma_tilde2 = tilde_func(sigma, 2*lambda)
      my_difference2 = round(abs(x) - abs(sigma_tilde2)/2, 10)
      condition2 = ((my_difference2 < 0) | ((my_difference2 == 0) & (sigma_tilde2 > 0)) | (x == 0))
      y = (1/lambda)*(1-sigma_tilde2/sigma) + (condition2)*(-1)^(sigma_tilde2 < 0)/sigma
      return(y)
    }
    return(sum_function)
  } else if(type == "linear") {
    sum_function = function(x) {
      x = tilde_func(x, lambda)
      sigma_tilde = tilde_func(sigma, lambda)
      my_difference = round(abs(x) - abs(sigma_tilde), 10)
      condition = ((my_difference < 0) | ((my_difference == 0) & (sigma_tilde > 0)) | (x == 0))
      y = (1/lambda)*(1-(sigma_tilde^2/sigma^2)) + (condition)*(abs(sigma_tilde)-abs(x))/sigma^2
      return(y)
    }
    return(sum_function)
  } else if(type == "exponential") {
    sum_function = function(x) {
      x = tilde_func(x, lambda)
      C = lambda/sigma
      if(any(x/sigma > 709)) {
         stop("overflow for the computation, sigma is too small or x is too large")
      }
      (1/sigma)*cosh(C-2*abs(x)/sigma)/sinh(C)
    }
    return(sum_function)
  } else if(type == "polynomial") {
    sum_function = function(x) {
      A = 2*pi*x/lambda
      B = 2*sigma/lambda
      (1/lambda)*sinh(B)/(cosh(B)-cos(A))
    }
    return(sum_function)
  } else if(type == "gaussian") {
    return(NULL)
  } else if(type == "sinc") {
    sum_function = function(x) {
      if(x %% lambda == 0) {
        # continuous in 0 in all cases
        return((1/lambda)*(2*floor(lambda/(2*sigma))+1))
      } else {
        return((1/(lambda*sin(pi*x/lambda)))*sin((2*floor(lambda/(2*sigma))+1)*pi*x/lambda))
      }
    }
    return(function(x) {sapply(x, sum_function)})
  } else if(type == "sinc2") {
    sum_function = function(x) {
      if(x %% lambda == 0) {
        N = floor(lambda/sigma)
        left = (1/lambda)*(2*floor(lambda/sigma)+1)
        right = -sigma*N*(N+1)/lambda^2
        return(left+right)
      } else {
        N = floor(lambda/sigma)
        left = (1/(lambda*sin(pi*x/lambda)))*sin((2*N+1)*pi*x/lambda)
        right = -sigma/(2*lambda^2)*((N+1)*cos(2*pi*N*x/lambda)-N*cos(2*pi*(N+1)*x/lambda)-1)/(sin(pi*x/lambda)^2)
        return(left+right)
      }
    }
    return(function(x) {sapply(x, sum_function)})
  } else {
    return(NULL)
  }
}

#' Sum of g (closed form)
#'
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc"
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @param lambda positive parameter conveying the periodicity
#' @returns sum of g taking a real `x` input and giving `Sg(x)` for this
#' `(type, sigma, lambda)` on its closed-form, or `NULL` if the closed-form
#' does not exist or is unknown
#' @export
Sg_closed_func = function(type, sigma, lambda) {
  if(is.null(type)) {
    return(NULL)
  } else if(type == "rectangular") {
    sum_function = function(x) {
      # the real Sg function is not defined in some points,
      # below it is set to 0 for those points, as for g_func
      # (if we would like to speak about distributions, we would have
      # either 0, +Dirac, -Dirac, or even Dirac' if Dirac and -Dirac are
      # at the same position)
      rep(0, length(x))
    }
    return(sum_function)
  } else if(type == "linear") {
    sum_function = function(x) {
      # the real Sg function is not defined in some points,
      # below it is set to 0 for those points, as for g_func
      x = tilde_func(x, lambda)
      sigma_tilde = tilde_func(sigma, lambda)
      my_difference = round(abs(x) - abs(sigma_tilde), 10)
      # we use the same condition as before for convenience
      condition = ((my_difference < 0) | ((my_difference == 0) & (sigma_tilde > 0)) | (x == 0))
      -condition*sign(x)/sigma^2
    }
    return(sum_function)
  } else if(type == "exponential") {
    sum_function = function(x) {
      x = tilde_func(x, lambda)
      C = lambda/sigma
      if(any(x/sigma > 709)) {
        stop("overflow for the computation, sigma is too small or x is too large")
      }
      # the real Sg function is not defined in some points,
      # below it is set to 0 for those points, as for g_func
      (-2*sign(x)/sigma^2)*sinh(C-2*abs(x)/sigma)/sinh(C)
    }
    return(sum_function)
  } else if(type == "polynomial") {
    sum_function = function(x) {
      A = 2*pi*x/lambda
      B = 2*sigma/lambda
      -(2*pi/lambda^2)*sin(A)*sinh(B)/(cosh(B)-cos(A))^2
    }
    return(sum_function)
  } else if(type == "gaussian") {
    return(NULL)
  } else if(type == "sinc") {
    sum_function = function(x) {
      A = lambda/2
      K = floor(A/sigma)
      L = K+1
      if(x %% lambda == 0) {
        return(0)
      } else {
        return(-pi/(lambda^2*(sin(pi*x/lambda))^2)*(L*sin(pi*x*K/A)-K*sin(pi*x*L/A)))
      }
    }
    return(function(x) {sapply(x, sum_function)})
  }  else if(type == "sinc2") {
    sum_function = function(x) {
      if(x %% lambda == 0) {
        return(0)
      } else {
        K = floor(lambda/sigma)
        L = K+1
        cos_diff = K*(cos(K*2*pi*x/lambda)-cos(L*2*pi*x/lambda))
        sin_diff = K*(sin(K*2*pi*x/lambda)-sin(L*2*pi*x/lambda))
        left = -pi/(lambda^2*(sin(pi*x/lambda))^2)*(L*sin(K*2*pi*x/lambda)-K*sin(L*2*pi*x/lambda))
        right = (pi*sigma/lambda^3)*(cos(pi*x/lambda)*(cos(K*2*pi*x/lambda)+cos_diff-1)+sin(pi*x/lambda)*L*sin_diff)/sin(pi*x/lambda)^3
        return(left+right)
      }
    }
    return(function(x) {sapply(x, sum_function)})
  } else {
    return(NULL)
  }
}
