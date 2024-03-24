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
#' @param q Input number or vector
#' @param maxiter Maximum number of iterations used. Note that the series
#' generally converge very quickly
#' @returns the numeric value phi(q)
#' @export
phi = function(q, maxiter = 3e4) {
  # function used to get explicit values for the Gaussian type
  # theta3(z, q)
  theta3(0, q = q, maxiter = maxiter)
}

#' theta3 function, similar to the method used in the elliptic package,
#' but ensuring that the maxiter doesn't stop after step 1
#' (e.g. elliptic::theta3(z = pi/4, q = exp(-pi/4), maxiter = 10000)
#' gives 1 whereas it should give 0.9135722, such as observed for
#' elliptic::theta3(z = pi/4+1e-15, q = exp(-pi/4), maxiter = 10000) ).
#' This is not sufficient because there can be some initial values close to
#' 0 at first (e.g. elliptic::theta3(z = pi/4, q = exp(-pi/16), maxiter = 10000)
#' gives 1 whereas it should give 0.1745516, such as observed for
#' elliptic::theta3(z = pi/4+1e-15, q = exp(-pi/16), maxiter = 10000) ).
#' We add below a comparison to the value before (may not solve the problem
#' in all the cases)
#'
#' @param z Primary complex argument
#' @param ignore Dummy argument to force the user to name the next
#' argument either `m` or `q`
#' @param m m as documented in `elliptic::theta3`
#' @param q q as documented in `elliptic::theta3`
#' @param give.n Boolean with default `FALSE` meaning to return the
#' function evaluation, and `TRUE` meaning to return a two element
#' list, with first element the function evaluation, and second element
#' the number of iterations used
#' @param maxiter Maximum number of iterations used. Note that the series
#' generally converge very quickly
#' @returns the numeric value theta3(z,q)
#' @export
theta3 = function(z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE,
          maxiter = 30) {
  if (!xor(is.null(m), is.null(q))) {
    stop("supply exactly one of m, q")
  }
  if (is.null(q)) {
    q <- elliptic::nome(m)
  }
  out <- 0
  out_previous <- -1
  for (n in 1:maxiter) {
    out.new <- out + q^(n^2) * cos(2 * z * n)
    if (elliptic::near.match(out, out.new) & n >= 5 & elliptic::near.match(out, out_previous)) {
      ans <- 1 + 2 * out
      if (give.n) {
        return(list(iterations = n, ans = ans))
      }
      else {
        return(ans)
      }
    }
    out_previous = out
    out <- out.new
  }
  stop("maximum iterations reached")
}

#' First derivative of the theta3 function, similar to the method used for the
#' derivative of theta1 found in the elliptic package. Additionally the same
#' correction is applied as for the theta3 function above
#'
#' @param z Primary complex argument
#' @param ignore Dummy argument to force the user to name the next
#' argument either `m` or `q`
#' @param m m as documented in `elliptic::theta3`
#' @param q q as documented in `elliptic::theta3`
#' @param give.n Boolean with default `FALSE` meaning to return the
#' function evaluation, and `TRUE` meaning to return a two element
#' list, with first element the function evaluation, and second element
#' the number of iterations used
#' @param maxiter Maximum number of iterations used. Note that the series
#' generally converge very quickly
#' @returns the numeric value theta3'(z,q)
#' @export
theta3dash = function (z, ignore = NULL, m = NULL, q = NULL, give.n = FALSE,
                       maxiter = 30) {
  if (!xor(is.null(m), is.null(q))) {
    stop("supply exactly one of m, q")
  }
  if (is.null(q)) {
    q <- elliptic::nome(m)
  }
  out <- 0
  out_previous <- -1
  for (n in 1:maxiter) {
    out.new <- out - q^(n^2) * (2 * n) * sin(2 * z * n)
    if (elliptic::near.match(out, out.new) & n >= 5 & elliptic::near.match(out, out_previous)) {
      ans <- 2 * out # the derivative removed the 1
      if (give.n) {
        return(list(iterations = n, ans = ans))
      }
      else {
        return(ans)
      }
    }
    out_previous = out
    out <- out.new
  }
  stop("maximum iterations reached")
}

#' Tilde function, to push any number x or sigma into [-lambda/2, lambda/2)
#'
#' @param x Input number (either x or sigma)
#' @param lambda lambda positive parameter corresponding to the periodicity
#' @returns x_tilde the number in [-lambda/2, lambda/2) such that
#' `x = x_tilde+i*lambda` (with i integer)
#' @export
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

#' Fractional part, with output in [0, 1)
#'
#' @param t Input number
#' @returns {t}_{+} the number in [0, 1)
#' @export
frac_plus = function(t) {
  t - floor(t)
}

#' Fractional part shifted, with output in [-1/2, 1/2)
#'
#' @param t Input number
#' @returns {t}_{-} the number in [-1/2, 1/2)
#' @export
frac_minus = function(t) {
  tilde_func(t, 1)
}

#' Find the maximum of an implicit function, with some assumptions
#'
#' @param F_tz function that input t,z and output a numeric value.
#' The function is such that, for each t, there is one and only one root z_root(t),
#' and, over t, the function z_root has a maximum z0 on the interval t_interval. There
#' should be a unique maximum on the interval
#' @param t_interval interval of t
#' @param z_interval interval of z
#' @param verbose whether to show plot at each step
#' @param precBits precision of the numeric representation
#' @param length.out interval length where to find the maximum
#' @returns the detected interval where the maximum z0 is, and the
#' corresponding time interval where the true t0 value is.
#' @export
find_maximum_implicit = function(F_tz,
                        t_interval = c(1, 3),
                        z_interval = c(1e-8, 0.5),
                        verbose = TRUE,
                        precBits=64,
                        length.out=11) {
  t_range = seq(from = mpfr(t_interval[1], precBits = precBits),
                to = mpfr(t_interval[2], precBits = precBits),
                length.out = length.out)
  stop_recursion = FALSE
  iteration = 1
  while(!stop_recursion) {
    z_range = list()
    for(k in 1:length(t_range)) {
      t = t_range[[k]]
      f = function(z) {
        F_tz(t, z)
      }
      z_range[[k]] = try(unirootR(f, z_interval, tol = .Machine$double.eps), silent=TRUE)$root
    }
    z_range = do.call(c, z_range)
    if(verbose) {
      # note: after some iterations, it is too small for 32 bits precision
      plot(sapply(t_range, asNumeric), sapply(z_range, asNumeric), main = iteration)
      iteration = iteration + 1
    }
    idx_neg = which(diff(z_range) < 0)
    idx_pos = which(diff(z_range) > 0)

    if(length(idx_neg) == 0 || length(idx_pos) == 0 || max(idx_pos) >= min(idx_neg)) {
      stop_recursion = TRUE
    } else { # increasing then decreasing
      idx = max(idx_pos)
      t_range_old = t_range
      z_range_old = z_range
      t_range = seq(from = t_range[[idx]],
                    to = t_range[[idx+2]],
                    length.out = length.out)
    }
  }

  # keep the interval
  l = list(t=range(t_range_old), z=range(z_range_old))
  return(l)
}

#' Convert an interval of [t1,t2] with corresponding [z1,z2] to a single
#' string truncating the common decimals
#'
#' @param l a list with elements t and z, each being a range
#' @returns a vector of string with t then z
#' @export
maximum_number_to_str = function(l) {
  t1 = trimws(capture.output(str(l$t[[1]], give.head=FALSE, digits.d=100)))
  t2 = trimws(capture.output(str(l$t[[2]], give.head=FALSE, digits.d=100)))
  z1 = trimws(capture.output(str(l$z[[1]], give.head=FALSE, digits.d=100)))
  z2 = trimws(capture.output(str(l$z[[2]], give.head=FALSE, digits.d=100)))

  max_z = max(which(sapply(1:min(nchar(z1), nchar(z2)), function(i){substr(z1, 1, i) == substr(z2, 1, i)})))
  max_t = max(which(sapply(1:min(nchar(t1), nchar(t2)), function(i){substr(t1, 1, i) == substr(t2, 1, i)})))

  final_z = substr(z1, 1, max_z)
  final_t = substr(t1, 1, max_t)
  return(c(final_t, final_z))
}
