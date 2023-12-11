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

save_video_from_png_folder = function(my_folder,
                                      stop_positions_frame_idx = c(),
                                      stop_positions_duration_s = 3,
                                      video_wo_stop_time_s = 10) {
  files = list.files(my_folder, "*.png", full.names = TRUE)
  N = length(files)
  framerate = floor(N/video_wo_stop_time_s)
  dur_each = rep(1, N)
  if(length(stop_positions_frame_idx) > 0) {
    dur_each[stop_positions_frame_idx] = stop_positions_duration_s*framerate
  }
  av::av_encode_video(rep(files, dur_each),
                      framerate = framerate,
                      output = paste0(my_folder, ".mp4"))
}
