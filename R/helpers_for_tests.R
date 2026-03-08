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

#' Test of a generic formula, only used for tests of Sf
#'
#' @param formula A function to test
#' @param type string either "linear", "exponential", "polynomial", "gaussian",
#' or "sinc", defining the function `f`
#' @param original_function ground truth function parametrized with `type`
#' @param boundary_step if NULL, the N_tested parameter is used. Otherwise,
#' create a grid with steps of `boundary_step`, to take e.g. 1, 1/2, 1/4, 1/8
#' to finer and finer granularity
#' @param N_tested if boundary_step is NULL, create random values for x, sigma,
#' sigma, and test the formula for `N_tested` repetitions
#' @param Sf_approx which approximation is chosen, use Sf_approx_direct_func for
#' linear and exponential since it converges quickly
#' @returns pass the test is the formula in valid in those cases
test_formula = function(formula, type = "linear", original_function = f_func,
                        boundary_step = NULL, N_tested = 1000,
                        Sf_approx = Sf_approx_direct_func) {
  if(is.null(boundary_step)) {
    N = N_tested
    set.seed(1)
    df = data.frame(x = runif(N, -3, 3),
                    lambda = abs(runif(N, -3, 3)),
                    sigma = abs(runif(N, -3, 3)))
  } else {
    step = boundary_step
    x = seq(from = -2, to = 2, by = step)
    lambda = seq(from = step, to = 2, by = step)
    sigma = seq(from = step, to = 2, by = step)
    df = expand.grid(x = x, lambda = lambda, sigma = sigma)
  }
  N = nrow(df)
  Sf_x_approx = rep(NA, N)
  Sf_x_formula = rep(NA, N)
  for(i in 1:N) {
    x = df$x[i]
    lambda = df$lambda[i]
    sigma = df$sigma[i]
    interval = seq(from=-10000L, to = 10000L, by = 1L) # can still be different from the closed formula for very small values
    Sf_x_approx[i] = Sf_approx(original_function(type, sigma), lambda, interval)(x)
    Sf_x_formula[i] = formula(x, sigma, lambda)

    if(N > 100) {
      print(i)
    }
    if(abs(Sf_x_approx[i] - Sf_x_formula[i]) > 1e-5) {
      print(paste(lambda, sigma, x))
    }
    expect_equal(Sf_x_approx[i], Sf_x_formula[i])
  }
}

#' Other test of a generic formula, only used for tests
#'
#' @param f predicted function, that should be a function of (x, lambda, sigma)
#' @param g ground truth function, that should be a function of (x, lambda, sigma)
#' @param step grid over x, lambda, sigma with steps `step`, to take e.g. 1, 1/2, 1/4, 1/8,
#' to finer and finer granularity
#' @param x_range range of the interval for x
#' @param lambda_max max of the interval for lambda
#' @param sigma_max max of the interval for sigma
#' @returns check the validity of f == g on the grid of values
check_equal_func = function(f, g, step = 1/2^3, x_range = c(-2, 2), lambda_max = 2, sigma_max = 2) {
  x = seq(from = x_range[1], to = x_range[2], by = step)
  lambda = seq(from = step, to = lambda_max, by = step)
  sigma = seq(from = step, to = sigma_max, by = step)
  df = expand.grid(x = x, lambda = lambda, sigma = sigma)
  N = nrow(df)
  for(i in 1:N) {
    if(N > 1000) {
      if(i %% 100 == 0) {
        print(paste0(i, "/", N))
      }
    }
    predicted = f(df$x[i], df$lambda[i], df$sigma[i])
    groundtruth = g(df$x[i], df$lambda[i], df$sigma[i])
    expect_equal(predicted, groundtruth)
  }
}

#' Other test of a generic formula with (t,z), only used for tests
#'
#' @param f predicted function, that should be a function of (t, z)
#' @param g ground truth function, that should be a function of (t, z)
#' @param step grid over t, z, and sigma with steps `step`, to take e.g. 1, 1/2, 1/4, 1/8,
#' to finer and finer granularity
#' @param z_range range of the interval for z
#' @param t_max max of the interval for lambda
#' @returns check the validity of f == g on the grid of values
#' @export
check_equal_tz_func = function(f, g, step = 1/2^3, z_range = c(-2, 2), t_max = 2, tol = NULL) {
  t = seq(from = step, to = t_max, by = step)
  z = seq(from = z_range[1], to = z_range[2], by = step)
  df = expand.grid(z = z, t = t)
  N = nrow(df)
  for(i in 1:N) {
    if(N > 1000) {
      if(i %% 100 == 0) {
        print(paste0(i, "/", N))
      }
    }
    predicted = f(df$t[i], df$z[i])
    groundtruth = g(df$t[i], df$z[i])
    if(is.null(tol)) {
      expect_equal(predicted, groundtruth)
    } else {
      expect_equal(predicted, groundtruth, tolerance = tol)
    }
  }
}

#' Convert a folder of png into an mp4 video
#'
#' @param my_folder folder containing images, ordered in the right order
#' @param stop_positions_frame_idx specific positions where the video should pause
#' @param stop_positions_duration_s duration of the pause in the stop positions, in second
#' @param video_wo_stop_time_s whole duration of the video in seconds
#' @returns nothing, save into the parent of `my_folder`, with the name `my_folder.mp4`
#' @export
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

#' Region with lobes for the sinc2 graph (used three times in test-dynamics.R)
#'
#' @param n_t number of t values in the [0,1] interval
#' @param n_z number of z values in the [-0.5, 0.5] interval
#' @returns the plot with the lobes for sinc2
#' @export
background_region_sinc2 = function(n_t=2000, n_z=2000, default = TRUE, t_min=0.003, t_max=1.0, z_min=-0.5, z_max=0.5) {
  # --- Compute Z_t(z) = 2t * sum_{k=1}^{floor(1/t)} (1 - k*t) * cos(2*k*pi*z) ---
  # Z_sincsq = function(t, z) {
  #   N = floor(1/t)
  #   if (N < 1) return(0)
  #   k = 1:N
  #   2 * t * sum((1 - k * t) * cos(2 * k * pi * z))
  # }

  # --- Build evaluation grid ---
  t_vals = seq(t_min, t_max, length.out = n_t)
  z_vals = seq(z_min, z_max, length.out = n_z)

  # Vectorized computation (vectorize over z for each t)
  Z_grid = matrix(0, nrow = n_t, ncol = n_z)
  for (i in seq_along(t_vals)) {
    tt = t_vals[i]
    N = floor(1/tt)
    if (N < 1) next
    k = 1:N
    weights = 1 - k * tt
    # Vectorized: outer product gives (N x n_z) matrix of cos values
    cos_mat = cos(2 * pi * outer(k, z_vals))
    Z_grid[i, ] = 2 * tt * colSums(weights * cos_mat)
  }

  # --- Build data frame for raster fill (+ zones in gray) ---
  df_grid = expand.grid(t = t_vals, z = z_vals)
  df_grid$positive = as.vector(Z_grid) > 0

  # --- Contour lines at Z = 0 ---
  # Use contourLines from base R
  cl = contourLines(t_vals, z_vals, Z_grid, levels = 0)
  df_contour = do.call(rbind, lapply(seq_along(cl), function(i) {
    data.frame(t = cl[[i]]$x, z = cl[[i]]$y, group = i)
  }))

  # --- Vertical band boundaries ---
  N_max_band = min(floor(1/min(t_vals)), 200)
  t_boundaries = 1 / (2:N_max_band)
  t_boundaries = t_boundaries[t_boundaries >= min(t_vals) & t_boundaries <= max(t_vals)]

  # --- Plot ---
  p = ggplot() +
    # Gray fill for positive zones
    geom_raster(
      data = df_grid[df_grid$positive, ],
      aes(x = t, y = z),
      fill = "#EEEEEE"
    ) +
    # Zero contour curves
    geom_path(
      data = df_contour,
      aes(x = t, y = z, group = group),
      color = "darkgray", linewidth = 0.5
    )

  if(default) {
    p = p +
      # Band boundary vertical lines
      #geom_vline(
      #  xintercept = t_boundaries,
      #  color = "darkgray", linewidth = 0.2, alpha = 0.3
      #) +
      # Labels
      annotate("text", x = 0.375, y = 0, label = "+",
               fontface = 2, size = 5, hjust = "center", vjust = "middle") +
      annotate("text", x = 0.375, y = 0.4, label = "\u2212",
               fontface = 2, size = 5, hjust = "center", vjust = "middle") +
      annotate("text", x = 0.375, y = -0.4, label = "\u2212",
               fontface = 2, size = 5, hjust = "center", vjust = "middle") +
      theme_bw() +
      xlab("t") +
      ylab("z") +
      scale_x_continuous(
        breaks = c(0, 1/2, 1),
        minor_breaks = c(1/4, 3/4),
        labels = c("0", "1/2", "1")
      ) +
      scale_y_continuous(
        breaks = c(-1/2, 0, 1/2),
        minor_breaks = c(-1/4, 1/4),
        labels = c("-1/2", "0", "1/2"),
        limits = c(-1/2, 1/2)
      ) +
      coord_fixed(xlim = c(0, 1))
  }
  p
}
