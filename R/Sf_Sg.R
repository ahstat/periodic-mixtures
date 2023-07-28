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
  } else if(type == "linear") {
    sum_function = function(x) {
      x = x %% lambda
      floorminus = floor((sigma - x)/lambda)
      floorplus = floor((sigma + x)/lambda)
      y = (1/sigma)*(
        floorplus*(1 - (-x + lambda*floorplus/2 + lambda/2)/sigma) +
          (1+floorminus)*(1 - (x + lambda*floorminus/2 + lambda*as.numeric(x>sigma))/sigma)
      )
      # x2 = x %% lambda
      #OK:
      # y = (1/sigma)*(1 + floor((sigma - x2)/lambda) + floor((sigma + x2)/lambda) -
      #                  (x2/sigma)*(1+floor((sigma - x2)/lambda) - floor((sigma + x2)/lambda)) -
      #                  (lambda/sigma)*(1/2)*((floor((sigma + x2)/lambda)*(floor((sigma + x2)/lambda)+1))+(pmax(0, floor((sigma - x2)/lambda))*(pmax(0, floor((sigma - x2)/lambda))+1))))

      # #New:
      # floorminus = floor((sigma - x2)/lambda)
      # floorplus = floor((sigma + x2)/lambda)
      # y = (1/sigma)*(
      #   1 +
      #     floorminus +
      #     floorplus -
      #     (x2/sigma)*(1+floorminus - floorplus) -
      #     (lambda/sigma)*
      #     (
      #       (floorplus*(floorplus+1))/2 +
      #         ((as.numeric(x2 > sigma) + floorminus)*(as.numeric(x2 > sigma) + floorminus+1))/2
      #     )
      # )

      # #New2:
      # floorminus = floor((sigma - x2)/lambda)
      # floorplus = floor((sigma + x2)/lambda)
      # y = (1/sigma)*(
      #   1 +
      #     floorminus +
      #     floorplus -
      #     (x2/sigma)*(1+floorminus - floorplus) -
      #     (lambda/sigma)*
      #     (
      #       (floorplus*(floorplus+1))/2 + floorminus*(floorminus+1)/2 +
      #         as.numeric(x2 > sigma) *(1+floorminus)
      #         #((as.numeric(x2 > sigma) + floorminus)*(as.numeric(x2 > sigma) + floorminus + 1))/2
      #     )
      # )

      # #New3:
      # floorminus = floor((sigma - x2)/lambda)
      # floorplus = floor((sigma + x2)/lambda)
      # y = (1/sigma)*(
      #   1 +
      #     floorminus +
      #     floorplus -
      #     (x2/sigma)*(1+floorminus - floorplus) -
      #     (lambda/sigma)*
      #     (
      #       (floorplus*(floorplus+1))/2 + floorminus*(floorminus+1)/2 +
      #         as.numeric(x2 > sigma) *(1+floorminus)
      #       #((as.numeric(x2 > sigma) + floorminus)*(as.numeric(x2 > sigma) + floorminus + 1))/2
      #     )
      # )

      # if(!all(pmax(0, floorminus) == as.numeric(x2 > sigma) + floorminus)) {
      #   stop("Not pmax ok")
      # }

      # #New4:
      # floorminus = floor((sigma - x2)/lambda)
      # floorplus = floor((sigma + x2)/lambda)
      # y = (1/sigma^2)*(
      #   (1+floorminus)*(sigma - x2 - lambda*floorminus/2 - lambda*as.numeric(x2>sigma)) +
      #     floorplus*(sigma + x2 - lambda*floorplus/2 - lambda/2)
      # )

      #New5:
      # floorminus = floor((sigma - x2)/lambda)
      # floorplus = floor((sigma + x2)/lambda)
      # y = (1/sigma)*(
      #   (1+floorminus)*(1 - (x2 + lambda*floorminus/2 + lambda*as.numeric(x2>sigma))/sigma) +
      #     floorplus*(1 - (-x2 + lambda*floorplus/2 + lambda/2)/sigma)
      # )
      return(y)
    }
    return(sum_function)
  } else if(type == "exponential") {
    sum_function = function(x) {
      x = x %% lambda
      # (1/(2*sigma)) * (exp(-x/sigma) + 2 * cosh(x/sigma) / (exp(lambda/sigma) - 1))
      if(any(x/sigma > 709)) {
        stop("Overflow for the computation, sigma is too small or x is too large")
      }
      (1/(2*sigma)) * (exp(-x/sigma) / (1-exp(-lambda/sigma)) - exp(x/sigma) / (1-exp(lambda/sigma)))
    }
    return(sum_function)
  } else if(type == "polynomial") {
    return(NULL)
  } else if(type == "gaussian") {
    return(NULL)
  } else if(type == "sinc") {
    # sum_function = function(x) {
    #   N = floor(lambda / (2*pi*sigma))
    #   if(N < 1) {
    #     sum_part = 0
    #   } else {
    #     sum_part = sum(cos(2*pi*(1:N)*x/lambda))
    #   }
    #   (1/lambda) + (2/lambda) * sum_part
    # }
    # return(function(x) {sapply(x, sum_function)})

    sum_function = function(x) {
      K = floor(lambda/(2*pi*sigma))
      y = pi*x/lambda
      term1 = cos((K+1)*y)
      term2 = sin(K*y)
      term3 = sin(y)
      sum_part = term1 * term2 / term3
      if(sum(is.nan(sum_part)) > 0) {
        sum_part[is.nan(sum_part)] = K # when x small to 0
      }
      if(sum(is.na(sum_part)) > 0) {
        sum_part[is.na(sum_part)] = K
      }
      (1/lambda) + (2/lambda) * sum_part
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
  } else if(type == "linear") {
    sum_function = function(x) {
      x = x %% lambda
      (-1/sigma^2) * (1 + floor((sigma-x)/lambda) - floor((sigma+x)/lambda))
    }
    return(sum_function)
  } else if(type == "exponential") {
    sum_function = function(x) {
      x = x %% lambda
      -(1/(2*sigma^2)) * (exp(-x/sigma) / (1-exp(-lambda/sigma)) + exp(x/sigma) / (1-exp(lambda/sigma)))
    }
    return(sum_function)
  } else if(type == "polynomial") {
    return(NULL)
  } else if(type == "gaussian") {
    return(NULL)
  } else if(type == "sinc") {
    # sum_function = function(x) {
    #   N = floor(lambda / (2*pi*sigma))
    #   if(N < 1) {
    #     sum_part = 0
    #   } else {
    #     sum_part = sum( (1:N) * sin(2*pi*(1:N)*x/lambda))
    #   }
    #   -((4*pi)/lambda^2) * sum_part
    # }
    sum_function = function(x) {
      K = floor(lambda/(2*pi*sigma))
      y = 2*pi*x/lambda
      sum_part = (1/4)*(1/(sin(y/2)^2))*((K+1)*sin(K*y) - K*sin((K+1)*y))
      if(sum(is.nan(sum_part)) > 0) {
        sum_part[is.nan(sum_part)] = 0
      }
      if(sum(is.na(sum_part)) > 0) {
        sum_part[is.na(sum_part)] = 0
      }
      -((4*pi)/lambda^2) * sum_part
    }
    return(function(x) {sapply(x, sum_function)})
  } else {
    return(NULL)
  }
}
