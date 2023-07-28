#' Function f
#'
#' Function $f_{\sigma}(x)$.
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc"
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `x` input and giving `f(x)` for this `sigma`
#' @export
f_func = function(type, sigma) {
  if(type == "linear") {
    f = function(x) {
      sigma^(-1) * pmax(1-abs(x)/sigma, 0)
    }
  } else if(type == "exponential") {
    f = function(x) {
      (1/(2*sigma)) * exp(-abs(x)/sigma)
    }
  } else if(type == "polynomial") {
    f = function(x) {
      # C = sigma^(-1)*sinc(n^(-1))/2 # n > 1, n real
      # n = 2: C = sigma^(-1)/pi
      # C / (1+abs(x/sigma)^n)
      # n = 2:
      sigma^(-1) * pi^(-1) / (1+(x/sigma)^2)
    }
  } else if(type == "gaussian") {
    f = function(x) {
      (2*pi)^(-1/2) * sigma^(-1) * exp(-x^2/(2*sigma^2))
    }
  } else if(type == "sinc") {
    f = function(x) {
      out = sin(x/sigma)/(pi*x)
      if(sum(is.nan(out)) > 0) {
        out[is.nan(out)] = 1/(sigma*pi)
      }
      return(out)
    }
  } else {
    stop("Unknown type")
  }
  return(f)
}

#' Function g, the derivative of f
#'
#' Derivative $g_{\sigma}(x)$ of the function $f_{\sigma}(x)$.
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f`
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `x` input and giving `g(x) := f'(x)` for this `sigma`
#' @export
g_func = function(type, sigma) {
  if(type == "linear") {
    g = function(x) {
      -sign(x) * sigma^(-2) * as.numeric(abs(x) <= sigma)
    }
  } else if(type == "exponential") {
    g = function(x) {
      -sign(x) * (1/2) * sigma^(-2) * exp(-abs(x)/sigma)
    }
  } else if(type == "polynomial") {
    g = function(x) {
      # C = -sigma^(-n-1)*n*sinc(n^(-1))/2 # n > 1, n real
      # n = 2: C = -2*sigma^(-3)/pi
      # C * (x*abs(x)^(n-2)) / (1+abs(x/sigma)^n)^2
      # n = 2:
      -2 * sigma^(-3) * pi^(-1) * x / (1+(x/sigma)^2)^2
    }
  } else if(type == "gaussian") {
    g = function(x) {
      -(2*pi)^(-1/2) * sigma^(-3) * exp(-x^2/(2*sigma^2)) * x
    }
  } else if(type == "sinc") {
    g = function(x) {
      out = (x*cos(x/sigma) - sigma * sin(x/sigma))/(sigma*pi*x^2)
      if(sum(is.nan(out)) > 0) {
        out[is.nan(out)] = 0
      }
      return(out)
    }
  } else {
    stop("Unknown type")
  }
  return(g)
}

#' Fourier transform of f
#'
#' Fourier transform \mathcal{F}f_{\sigma}(\xi).
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f`
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `\xi` input and giving `\mathcal{F}f_{\sigma}(\xi)`
#' the Fourier transform of `f` in `\xi`
#' @export
Ff_func = function(type, sigma) {
  if(type == "linear") {
    Ff = function(xi) {
      (sinc(sigma * xi))^2
    }
  } else if(type == "exponential") {
    Ff = function(xi) {
      1/(1 + (2*pi*sigma*xi)^2)
    }
  } else if(type == "polynomial") {
    Ff = function(xi) {
      exp(-2*pi*sigma*abs(xi))
    }
  } else if(type == "gaussian") {
    Ff = function(xi) {
      exp(-(2*pi*sigma*xi)^2/2)
    }
  } else if(type == "sinc") {
    Ff = function(xi) {
      as.numeric(abs(xi) <= (1/(2*pi*sigma)))
    }
  } else {
    stop("Unknown type")
  }
  return(Ff)
}

#' Fourier transform of g
#'
#' Fourier transform \mathcal{F}g_{\sigma}(\xi) of g, the derivative of f.
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f`
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `\xi` input and giving `\mathcal{F}g_{\sigma}(\xi)`
#' the Fourier transform of `f` in `\xi`
#' @export
Fg_func = function(type, sigma) {
  Ff = Ff_func(type, sigma)
  Fg = function(xi) {
    2 * pi * 1i * xi * Ff(xi)
  }
  return(Fg)
}

#' sigma such as the Fourier transform sums to one
#'
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f`
#' @returns function the positive number `\sigma` such that `\mathcal{F}f_{\sigma}`
#' sums to one
#' @export
sigma_such_as_Fourier_tranform_sums_to_one_func = function(type) {
  if(type == "linear") {
    sigma = 1
  } else if(type == "exponential") {
    sigma = 1/2
  } else if(type == "polynomial") {
    sigma = 1/pi
  } else if(type == "gaussian") {
    sigma = 1/sqrt(2*pi)
  } else if(type == "sinc") {
    sigma = 1/pi
  } else {
    stop("Unknown type")
  }
  return(sigma)
}
