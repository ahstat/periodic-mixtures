test_that("f sums to 1 for all types", {
  test_f_sums_to_one = function(type, eps = 1e-5, N = 400, length.out = 1e5, verbose = TRUE) {
    sigmas = sigmas_and_lambdas_func(N)$sigmas
    x0 = seq(from = -1e8, to = 1e8, length.out = length.out) # for very large variance functions
    x1 = seq(from = -10000, to = 10000, length.out = length.out) # for large variance functions
    x2 = seq(from = -10, to = 10, length.out = length.out) # for small variance functions
    x3 = seq(from = -1, to = 1, length.out = length.out) # for very small variance functions
    vals = rep(NA, N)
    for(k in 1:N) {
      if(k %% 10 == 0) {
        if(verbose) {
          print(k)
        }
      }
      sigma = sigmas[k]
      f = f_func(type, sigma)

      # deviation of the integral on a different interval (also different granularities)
      val0 = abs(sum(f(x0))*(x0[2]-x0[1]) - 1)
      val1 = abs(sum(f(x1))*(x1[2]-x1[1]) - 1)
      val2 = abs(sum(f(x2))*(x2[2]-x2[1]) - 1)
      val3 = abs(sum(f(x3))*(x3[2]-x3[1]) - 1)

      # Take the best
      vals[k] = min(val0, val1, val2, val3)
    }
    if(verbose) {
      print(max(vals))
    }
    expect_true(all(vals <= eps)) # the integral is very close to 1 for all
  }

  ## Quick checks (large epsilon)
  expect_error(test_f_sums_to_one("rectangular", eps = 6e-5, N = 10, verbose = FALSE), NA)
  expect_error(test_f_sums_to_one("linear", eps = 1e-5, N = 10, verbose = FALSE), NA)
  expect_error(test_f_sums_to_one("exponential", eps = 9e-4, N = 10, verbose = FALSE), NA)
  expect_error(test_f_sums_to_one("polynomial", eps = 2e-2, N = 10, verbose = FALSE), NA)
  expect_error(test_f_sums_to_one("gaussian", eps = 4e-12, N = 10, verbose = FALSE), NA)
  expect_error(test_f_sums_to_one("sinc", eps = 2e-3, N = 10, verbose = FALSE), NA)
  expect_error(test_f_sums_to_one("sinc2", eps = 1e-3, N = 10, verbose = FALSE), NA)

  # ## Longer checks (smaller epsilon)
  # expect_error(test_f_sums_to_one("rectangular", eps = 5e-5, length.out = 1e6, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("linear", eps = 2e-7, length.out = 1e6, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("exponential", eps = 2e-5, length.out = 1e6, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("polynomial", eps = 2e-3, length.out = 1e6, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("gaussian", eps = 4e-11, length.out = 1e6, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("sinc", eps = 7e-4, length.out = 1e6, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("sinc2", eps = 2e-4, length.out = 1e6, verbose = TRUE), NA)
  #
  # ## Very long checks (for polynomial, slow convergence...)
  # expect_error(test_f_sums_to_one("exponential", eps = 4e-7, length.out = 1e7, verbose = TRUE), NA)
  # expect_error(test_f_sums_to_one("polynomial", eps = 7e-4, length.out = 1e7, verbose = TRUE), NA)
})

test_that("g is the derivative of f for all types", {
  test_g_is_derivative_of_f = function(type, eps = 1e-5, N = 400, length.out = 1e5, verbose = TRUE) {
    sigmas = sigmas_and_lambdas_func(N)$sigmas
    x1 = seq(from = -1000, to = 1000, length.out = length.out) # for large variance functions
    x2 = seq(from = -10, to = 10, length.out = length.out) # for small variance functions
    vals = rep(NA, N)
    for(k in 1:N) {
      if(k %% 10 == 0) {
        if(verbose) {
          print(k)
        }
      }
      sigma = sigmas[k]
      f = f_func(type, sigma)
      g = g_func(type, sigma)

      # comparing g(x) with f(x)
      val1 = abs(diff(f(x1))/(x1[2]-x1[1]) - g(x1)[-1])
      # remove possible wrong values (cf abrupt slope change)
      val1 = sort(val1, decreasing = TRUE)[-c(1:10)]
      val1 = mean(val1)
      val2 = abs(diff(f(x2))/(x2[2]-x2[1]) - g(x2)[-1])
      # remove possible wrong values (cf abrupt slope change)
      val2 = sort(val2, decreasing = TRUE)[-c(1:10)]
      val2 = mean(val2)

      # Take the best
      vals[k] = min(val1, val2)
    }
    if(verbose) {
      print(max(vals))
    }
    expect_true(all(vals <= eps)) # f' computed from scratch is very close to g
  }

  expect_error(test_g_is_derivative_of_f("rectangular", eps = 3e-15, N = 10, verbose = FALSE), NA)
  expect_error(test_g_is_derivative_of_f("linear", eps = 3e-15, N = 10, verbose = FALSE), NA)
  expect_error(test_g_is_derivative_of_f("exponential", eps = 3e-4, N = 10, verbose = FALSE), NA)
  expect_error(test_g_is_derivative_of_f("polynomial", eps = 2e-4, N = 10, verbose = FALSE), NA)
  expect_error(test_g_is_derivative_of_f("gaussian", eps = 3e-4, N = 10, verbose = FALSE), NA)
  expect_error(test_g_is_derivative_of_f("sinc", eps = 3e-2, N = 10, verbose = FALSE), NA)
  expect_error(test_g_is_derivative_of_f("sinc2", eps = 2e-3, N = 10, verbose = FALSE), NA)

  # expect_error(test_g_is_derivative_of_f("linear", eps = 3e-15, verbose = TRUE), NA)
  # expect_error(test_g_is_derivative_of_f("exponential", eps = 3e-4, verbose = TRUE), NA)
  # expect_error(test_g_is_derivative_of_f("polynomial", eps = 3e-4, verbose = TRUE), NA)
  # expect_error(test_g_is_derivative_of_f("gaussian", eps = 4e-4, verbose = TRUE), NA)
})

test_that("sigma as defined by `sigma_such_as_Fourier_tranform_sums_to_one_func` makes `Ff_{sigma}` sum to one", {
  test_sigma_in_table_gives_Ff_sums_to_one = function(type, eps = 1e-5, length.out = 1e5, verbose = TRUE) {
    sigma = sigma_such_as_Fourier_tranform_sums_to_one_func(type)

    x0 = seq(from = -1e8, to = 1e8, length.out = length.out) # for very large variance functions
    x1 = seq(from = -10000, to = 10000, length.out = length.out) # for large variance functions
    x2 = seq(from = -10, to = 10, length.out = length.out) # for small variance functions

    Ff = Ff_func(type, sigma)

    # deviation of the integral on a different interval (also different granularities)
    val0 = abs(sum(Ff(x0))*(x0[2]-x0[1]) - 1)
    val1 = abs(sum(Ff(x1))*(x1[2]-x1[1]) - 1)
    val2 = abs(sum(Ff(x2))*(x2[2]-x2[1]) - 1)

    # Take the best
    val = min(val0, val1, val2)

    if(verbose) {
      print(val)
    }
    expect_true(val <= eps) # the integral is very close to 1
  }

  ## Quick checks (large epsilon)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("rectangular", eps = 2e-5, verbose = FALSE), NA)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("linear", eps = 2e-5, verbose = FALSE), NA)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("exponential", eps = 2e-4, verbose = FALSE), NA)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("polynomial", eps = 9e-9, verbose = FALSE), NA)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("gaussian", eps = 3e-12, verbose = FALSE), NA)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("sinc", eps = 2e-5, verbose = FALSE), NA)
  expect_error(test_sigma_in_table_gives_Ff_sums_to_one("sinc2", eps = 2e-10, verbose = FALSE), NA)
})

test_that("plots for f, g, Ff, Fg output correctly", {
  recompute = FALSE
  output_folder_plots0 = "~/Documents/GitHub/ahstat.github.io/images"
  # output_folder_plots0 = "~/Github/ahstat.github.io/images"
  if(!dir.exists(output_folder_plots0)) {
    # pass this sequence, since there is no corresponding folder
    expect_true(0 == 0)
  } else {
    # library(periodicmixtures)
    library(ggplot2)
    library(ggh4x)
    output_folder_plots = file.path(output_folder_plots0, "2023-6-11-Periodic-mixtures/plot1")
    dir.create(output_folder_plots, showWarnings = FALSE)
    subfolders = list.files(output_folder_plots)
    if(identical(subfolders, c("f", "ℱf", "ℱg", "g")) & !recompute) {
      # pass this sequence, since the folder already exists
      # delete the folder if you want to recompute it
      expect_true(0 == 0)
    } else {
      # there are missing folders, recompute all of them
      types = c("rectangular", "linear", "exponential", "polynomial", "gaussian", "sinc", "sinc2")
      x = seq(from = -3, to = 3, length.out = 30001)
      p = list() # output plots
      my_functions = c("f" = f_func,
                       "g" = g_func,
                       "ℱf" = Ff_func,
                       "ℱg" = Fg_func)
      ylims = list(c(0,1), c(-3,3), c(0,1), c(-3,3)) # ylims, respective with `my_functions` elements
      freq_variable = "ξ"
      # compute each function for each type
      for(k in 1:length(my_functions)) {
        name = names(my_functions)[k]
        # print(name)
        my_function = my_functions[[k]]
        my_ylim = ylims[[k]]
        p[[name]] = list()
        for(type in types) {
          sigma = sigma_such_as_Fourier_tranform_sums_to_one_func(type)
          f = my_function(type, sigma)
          variable = ifelse(grepl("ℱ", name), freq_variable, "x")
          ylab_name = bquote(.(name)[σ](.(variable)))
          xlab_name = variable
          y = f(x)
          if(name == "ℱg" ) {
            if(!all(Re(f(x)) == 0)) {
              stop("there is a real part")
            }
            ylab_name = bquote(Im(.(name)[σ](.(variable)))) # paste0("Im(", ylab_name, ")")
            y = Im(y)
          }
          data = data.frame(x=x,y=y)
          p[[name]][[type]] = ggplot(data = data, aes(x=x,y=y)) +
            geom_line() +
            theme_bw() +
            ylab(ylab_name) +
            xlab(xlab_name) +
            coord_cartesian(ylim = my_ylim)
        }
      }
      # plot each function for each type
      for(name1 in names(p)) {
        for(name2 in names(p[[name1]])) {
          dir.create(file.path(output_folder_plots, name1), recursive = TRUE, showWarnings = FALSE)
          ggsave(filename = file.path(output_folder_plots, name1, paste0(name2, ".png")),
                 plot = p[[name1]][[name2]]+force_panelsizes(rows = unit(2, "cm"),cols = unit(2, "cm")),
                 width = 120, height = 120, units = "px", dpi = 85) # larger dpi = larger font
        }
      }
      expect_true(0 == 0)
    }
  }
})
