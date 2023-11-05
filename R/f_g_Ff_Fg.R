#' Function f
#'
#' Function $f_{\sigma}(x)$.
#' @param type string either "rectangular", "linear", "exponential", "polynomial",
#' "gaussian", "sinc", or "sinc2"
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `x` input and giving `f(x)` for this `sigma`
#' @export
f_func = function(type, sigma) {
  if(type == "rectangular") {
    f = function(x) {
      sigma^(-1) * as.numeric(abs(x) <= (sigma/2))
    }
  } else if(type == "linear") {
    f = function(x) {
      sigma^(-1) * pmax(1-abs(x)/sigma, 0)
    }
  } else if(type == "exponential") {
    f = function(x) {
      sigma^(-1) * exp(-2*abs(x)/sigma)
    }
  } else if(type == "polynomial") {
    f = function(x) {
      sigma^(-1) / (1+(pi*x/sigma)^2)
    }
  } else if(type == "gaussian") {
    f = function(x) {
      sigma^(-1) * exp(-pi*x^2/sigma^2)
    }
  } else if(type == "sinc") {
    f = function(x) {
      sigma^(-1) * sinc(x/sigma)
    }
  } else if(type == "sinc2") {
    f = function(x) {
      sigma^(-1) * sinc(x/sigma)^2
    }
  }
  else {
    stop("Unknown type")
  }
  return(f)
}

#' Function g, the derivative of f
#'
#' Derivative $g_{\sigma}(x)$ of the function $f_{\sigma}(x)$.
#' @param type string either "rectangular", "linear", "exponential", "polynomial",
#' "gaussian", "sinc", or "sinc2", defining the function `f`
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `x` input and giving `g(x) := f'(x)` for this `sigma`
#' @export
g_func = function(type, sigma) {
  if(type == "rectangular") {
    g = function(x) {
      # not defined in some points, below it is set to plus or minus 1 for those points
      out = rep(0, length(x))
      positive_peak = which(x == -sigma/2)
      negative_peak = which(x == sigma/2)
      if(length(positive_peak) > 0) {
        out[positive_peak] = 1
      }
      if(length(negative_peak) > 0) {
        out[negative_peak] = -1
      }
      return(out)
    }
  } else if(type == "linear") {
    g = function(x) {
      # not defined in some points, below it is set to 0 for those points
      -sign(x) * sigma^(-2) * as.numeric(abs(x) <= sigma)
    }
  } else if(type == "exponential") {
    g = function(x) {
      -sign(x) * 2 * sigma^(-2) * exp(-2*abs(x)/sigma)
    }
  } else if(type == "polynomial") {
    g = function(x) {
      -2 * pi^2 * sigma^(-3) * x / (1+(pi*x/sigma)^2)^2
    }
  } else if(type == "gaussian") {
    g = function(x) {
      -2 * pi * sigma^(-3) * x * exp(-pi*x^2/sigma^2)
    }
  } else if(type == "sinc") {
    g = function(x) {
      out = (1/(sigma*x))*cos(pi*x/sigma) - (1/(pi*x^2)) * sin(pi*x/sigma)
      if(sum(is.nan(out)) > 0) {
        out[is.nan(out)] = 0
      }
      return(out)
    }
  } else if(type == "sinc2") {
    g = function(x) {
      out = sin(2*pi*x/sigma)/(pi*x^2) - 2*sigma*sin(pi*x/sigma)^2/(pi^2*x^3)
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
#' @param type string either "rectangular", "linear", "exponential", "polynomial",
#' "gaussian", "sinc", or "sinc2", defining the function `f`
#' @param sigma positive parameter conveying variance information, such as the
#' sd for the gaussian type
#' @returns function taking a real `\xi` input and giving `\mathcal{F}f_{\sigma}(\xi)`
#' the Fourier transform of `f` in `\xi`
#' @export
Ff_func = function(type, sigma) {
  if(type == "rectangular") {
    Ff = function(xi) {
      sinc(sigma * xi)
    }
  } else if(type == "linear") {
    Ff = function(xi) {
      (sinc(sigma * xi))^2
    }
  } else if(type == "exponential") {
    Ff = function(xi) {
      1/(1 + (pi*sigma*xi)^2)
    }
  } else if(type == "polynomial") {
    Ff = function(xi) {
      exp(-2*sigma*abs(xi))
    }
  } else if(type == "gaussian") {
    Ff = function(xi) {
      exp(-pi*sigma^2*xi^2)
    }
  } else if(type == "sinc") {
    Ff = function(xi) {
      as.numeric(abs(xi) <= (1/(2*sigma)))
    }
  } else if(type == "sinc2") {
    Ff = function(xi) {
      pmax(1-abs(xi)*sigma, 0)
    }
  } else {
    stop("Unknown type")
  }
  return(Ff)
}

#' Fourier transform of g
#'
#' Fourier transform \mathcal{F}g_{\sigma}(\xi) of g, the derivative of f.
#' @param type string either "rectangular", "linear", "exponential", "polynomial",
#' "gaussian", "sinc", or "sinc2", defining the function `f`
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
#' (this has been aligned so that the sigma to select is always 1)
#'
#' @param type string either "rectangular", "linear", "exponential", "polynomial",
#' "gaussian", "sinc", or "sinc2", defining the function `f`
#' @returns function the positive number `\sigma` such that `\mathcal{F}f_{\sigma}`
#' sums to one
#' @export
sigma_such_as_Fourier_tranform_sums_to_one_func = function(type) {
  if(type == "rectangular") {
    sigma = 1
  } else if(type == "linear") {
    sigma = 1
  } else if(type == "exponential") {
    sigma = 1
  } else if(type == "polynomial") {
    sigma = 1
  } else if(type == "gaussian") {
    sigma = 1
  } else if(type == "sinc") {
    sigma = 1
  } else if(type == "sinc2") {
    sigma = 1
  } else {
    stop("Unknown type")
  }
  return(sigma)
}
