test_that("Sf_approx_direct and Sf_approx_Fourier are almost equal for large interval K, for all types", {
  test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same = function(type, eps = 1e-5, N = 10, length.out = 1e5, K = 100L, verbose = TRUE) {
    sigmas_and_lambdas = sigmas_and_lambdas_func(N)
    sigmas = sigmas_and_lambdas$sigmas
    lambdas = sigmas_and_lambdas$lambdas
    vals = rep(NA, N)
    for(k in 1:N) {
      if(k %% 1 == 0) {
        if(verbose) {
          print(k)
        }
      }
      sigma = sigmas[k]
      lambda = lambdas[k]
      interval = seq(from = -K, to = K, by = 1L)

      elements = build_f_g_Ff_Fg(type, sigma)
      f = elements$f
      Ff = elements$Ff

      Sf_approx_direct = Sf_approx_direct_func(f, lambda, interval)
      Sf_approx_Fourier = Sf_approx_Fourier_func(Ff, lambda, interval)

      # Interval depending of lambda the periodicity
      x = seq(from = -2*lambda, to = 1*lambda, length.out = length.out)
      if(verbose) {
        plot(x, Sf_approx_direct(x), type = "l", main = k)
        lines(x, Sf_approx_Fourier(x), col = "red")
      }
      val = mean(abs(Sf_approx_direct(x) - Sf_approx_Fourier(x)))

      vals[k] = val
    }
    if(verbose) {
      print(max(vals))
    }
    expect_true(all(vals <= eps)) # difference between direct sum and Fourier sum is small
    # if true, it's ok, direct sum and Fourier sum give (almost) the same
  }

  ## Quick checks
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("rectangular", eps = 4e-4, K = 1000L, N=1, length.out=1e3, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("linear", eps = 7e-3, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("exponential", eps = 2e-3, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("gaussian", eps = 2e-15, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("polynomial", eps = 2e-3, N=1, length.out=1e2, K = 1000L, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc", eps = 2e-2, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc", eps = 8e-4, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc2", eps = 3e-6, K = 1000L, N=1, length.out=1e3, verbose = FALSE), NA)
  # ## Long checks (200 seconds for 6*2 tests)
  # expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("linear", eps = 7e-3, verbose = TRUE), NA)
  # expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("exponential", eps = 2e-3, verbose = TRUE), NA)
  # expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("gaussian", eps = 2e-15, verbose = TRUE), NA)
  # # for the three following, sometimes looks different (with verbose=TRUE plots)
  # # because almost constant function (but it should work)
  # expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("polynomial", eps = 2e-3, K = 1000L, verbose = TRUE), NA)
  # expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc", eps = 2e-2, verbose = TRUE), NA)
  # expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc", eps = 8e-4, K = 1000L, verbose = TRUE), NA)
})

test_that("Sg_approx_direct and Sg_approx_Fourier are almost equal for large interval K, for all types", {
  test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same = function(type, eps = 1e-5, N = 10, length.out = 1e5, K = 100L, verbose = TRUE) {
    sigmas_and_lambdas = sigmas_and_lambdas_func(N)
    sigmas = sigmas_and_lambdas$sigmas
    lambdas = sigmas_and_lambdas$lambdas
    vals = rep(NA, N)
    for(k in 1:N) {
      if(k %% 1 == 0) {
        if(verbose) {
          print(k)
        }
      }
      sigma = sigmas[k]
      lambda = lambdas[k]
      interval = seq(from = -K, to = K, by = 1L)

      elements = build_f_g_Ff_Fg(type, sigma)
      g = elements$g
      Ff = elements$Ff

      Sg_approx_direct = Sg_approx_direct_func(g, lambda, interval) # g needed
      Sg_approx_Fourier = Sg_approx_Fourier_func(Ff, lambda, interval) # Ff needed

      # Interval depending of lambda the periodicity
      x = seq(from = -2*lambda, to = 1*lambda, length.out = length.out)
      if(verbose) {
        plot(x, round(Sg_approx_direct(x), 10), type = "l", main = k)
        lines(x, round(Sg_approx_Fourier(x), 10), col = "red")
      }
      val = mean(abs(Sg_approx_direct(x) - Sg_approx_Fourier(x)))

      vals[k] = val
    }
    if(verbose) {
      print(max(vals))
    }
    expect_true(all(vals <= eps)) # difference between direct sum and Fourier sum is small
    # if true, it's ok, direct sum and Fourier sum give (almost) the same
  }

  ## Quick checks
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("rectangular", eps = 2, N=1, length.out=1e2,verbose = FALSE), NA) # not giving good results
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("linear", eps = 2, N=1, length.out=1e2, verbose = FALSE), NA) # large Gibbs effect
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("exponential", eps = 0.6, N=1, length.out=1e2, verbose = FALSE), NA) # large Gibbs effect
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("gaussian", eps = 2e-14, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("polynomial", eps = 2e-9, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("sinc", eps = 6e-4, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("sinc2", eps = 6e-9, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)

  # ## Long checks (240 seconds for 5*2 tests)
  # expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("linear", eps = 2, verbose = FALSE), NA) # large Gibbs effect
  # expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("exponential", eps = 0.6, verbose = FALSE), NA) # large Gibbs effect
  # expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("gaussian", eps = 2e-14, verbose = FALSE), NA)
  # expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("polynomial", eps = 2e-9, K = 1000L, verbose = FALSE), NA)
  # expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("sinc", eps = 6e-4, K = 1000L, verbose = FALSE), NA)
})

test_that("Sf_approx (either Fourier or direct depending on cases) and Sf_closed are almost equal for large interval K, for all types for which the closed-form is known", {
  test_Sf_approx_and_Sf_closed_gives_same = function(type, eps = 1e-5, N = 100, length.out = 1e5, K = 100L, verbose = TRUE, approx_is_Fourier = FALSE) {
    sigmas_and_lambdas = sigmas_and_lambdas_func(N)
    sigmas = sigmas_and_lambdas$sigmas
    lambdas = sigmas_and_lambdas$lambdas
    vals = rep(NA, N)
    for(k in 1:N) {
      if(k %% 1 == 0) {
        if(verbose) {
          print(k)
        }
      }
      sigma = sigmas[k]
      lambda = lambdas[k]
      interval = seq(from = -K, to = K, by = 1L)

      elements = build_f_g_Ff_Fg(type, sigma)

      if(approx_is_Fourier) {
        Ff = elements$Ff
        Sf_approx = Sf_approx_Fourier_func(Ff, lambda, interval)
      } else {
        f = elements$f
        Sf_approx = Sf_approx_direct_func(f, lambda, interval)
      }

      Sf_closed = Sf_closed_func(type, sigma, lambda)

      # Interval depending of lambda the periodicity
      x = seq(from = -2*lambda, to = 1*lambda, length.out = length.out)
      if(verbose) {
        plot(x, Sf_approx(x), type = "l", main = k)
        lines(x, Sf_closed(x), col = "red")
      }
      val = mean(abs(Sf_approx(x) - Sf_closed(x)))
      vals[k] = val
    }
    if(verbose) {
      print(max(vals))
    }
    expect_true(all(vals <= eps)) # difference between direct sum and closed form sum is small
    # if true, it's ok, direct sum and closed form sum give the same
  }

  ## Quick checks
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("rectangular", eps = 1e-16, K = 1000L, approx_is_Fourier = FALSE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("linear", eps = 4e-16, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("exponential", eps = 9e-11, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("polynomial", eps = 6e-5, N=1, length.out=1e2, verbose = FALSE), NA) # checked with N=10 and a larger K
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("gaussian", eps = 7e-15, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 1e-8, approx_is_Fourier = TRUE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 2e-1, K = 1000L, approx_is_Fourier = FALSE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc2", eps = 1e-8, approx_is_Fourier = TRUE, N=1, length.out=1e2, verbose = FALSE), NA)

  ## Long checks (~1800 seconds for 4*2 tests)
  # expect_error(test_Sf_approx_and_Sf_closed_gives_same("linear", eps = 4e-16, verbose = TRUE), NA)
  # expect_error(test_Sf_approx_and_Sf_closed_gives_same("exponential", eps = 9e-11, verbose = TRUE), NA)
  # # For exponential, k=10, (sigma=5.14, lambda=0.39) we have 1.409717e-15 difference with K=1000L but only 0.001362376 with K=100L
  # # For exponential, k=86, (sigma=8.99, lambda=0.22) we have 8.230421e-11 difference with K=1000L but only 0.3758261 with K=100L
  # expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 3e-15, approx_is_Fourier = TRUE, verbose = TRUE), NA)
  # expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 2e-1, K = 1000L, approx_is_Fourier = FALSE, verbose = TRUE), NA)
  # # expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 6e-3, K = 10000L, approx_is_Fourier = FALSE, verbose = TRUE), NA) # long
  # there is no known closed formulas for polynomial or gaussian
})

test_that("Sg_approx (either Fourier or direct depending on cases) and Sg_closed are almost equal for large interval K, for all types for which the closed-form is known", {
  test_Sg_approx_and_Sg_closed_gives_same = function(type, eps = 1e-5, N = 100, length.out = 1e4, K = 1000L, verbose = TRUE, approx_is_Fourier = FALSE) {
    sigmas_and_lambdas = sigmas_and_lambdas_func(N)
    sigmas = sigmas_and_lambdas$sigmas
    lambdas = sigmas_and_lambdas$lambdas
    vals = rep(NA, N)
    for(k in 1:N) {
      if(k %% 1 == 0) {
        if(verbose) {
          print(k)
        }
      }
      sigma = sigmas[k]
      lambda = lambdas[k]
      interval = seq(from = -K, to = K, by = 1L)

      elements = build_f_g_Ff_Fg(type, sigma)

      if(approx_is_Fourier) {
        Ff = elements$Ff
        Sg_approx = Sg_approx_Fourier_func(Ff, lambda, interval)
      } else {
        g = elements$g
        Sg_approx = Sg_approx_direct_func(g, lambda, interval)
      }

      Sg_closed = Sg_closed_func(type, sigma, lambda)

      # Interval depending of lambda the periodicity
      x = seq(from = -2*lambda, to = 1*lambda, length.out = length.out)
      if(verbose) {
        plot(x, Sg_approx(x), type = "l", main = k)
        lines(x, Sg_closed(x), col = "red")
      }
      val = abs(Sg_approx(x) - Sg_closed(x))
      val = sort(val, decreasing = TRUE)[-c(1:10)] # remove possible wrong values (cf abrupt slope change)
      val = mean(val)
      vals[k] = val
    }
    if(verbose) {
      print(max(vals))
    }
    expect_true(all(vals <= eps)) # difference between direct sum and Fourier sum is small
    # if true, it's ok, direct sum and Fourier sum give (almost) the same
  }

  ## Quick checks
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("rectangular", eps = 3e-16, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("linear", eps = 3e-16, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("exponential", eps = 2e-13, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("polynomial", eps = 2e-12, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("gaussian", eps = 2e-13, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc", eps = 4e-13, approx_is_Fourier = TRUE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc", eps = 8e-4, K = 1000L, approx_is_Fourier = FALSE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc2", eps = 4e-13, approx_is_Fourier = TRUE, N=1, length.out=1e2, verbose = FALSE), NA)

  # ## Long checks (437 seconds for 4*2 tests)
  # expect_error(test_Sg_approx_and_Sg_closed_gives_same("linear", eps = 3e-16, verbose = TRUE), NA)
  # expect_error(test_Sg_approx_and_Sg_closed_gives_same("exponential", eps = 2e-13, verbose = TRUE), NA)
  # expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc", eps = 4e-13, approx_is_Fourier = TRUE, verbose = TRUE), NA)
  # expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc", eps = 8e-4, K = 1000L, approx_is_Fourier = FALSE, verbose = TRUE), NA)
  # # there is no known closed formulas for polynomial or gaussian
})

test_that("videos and plots for Sf output correctly", {
  recompute = FALSE
  output_folder_plots0 = "~/Documents/GitHub/ahstat.github.io/images"
  if(!dir.exists(output_folder_plots0)) {
    # pass this sequence, since there is no corresponding folder
    expect_true(0 == 0)
  } else {
    # library(periodicmixtures)
    library(ggplot2)
    library(ggh4x)
    output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot2"
    dir.create(output_folder_plots, showWarnings = FALSE)
    subfolders = list.files(output_folder_plots)
    if(("Sf" %in% subfolders) & !recompute) {
      # pass this sequence, since the folder already exists
      # delete the folder if you want to recompute it
      expect_true(0 == 0)
    } else {
      # there are missing folders, recompute all of them

      ## Step 1. Plot all the frames for the video for Sf ----
      plot_Sf = function(type, sigma, lambda,
                         to_plot = c("closed", "direct", "Fourier"),
                         length_x = 3001, K = 3000L, add_lambda = TRUE) {
        x = seq(from = -lambda, to = lambda, length.out = length_x)
        x = x[-c(1, length(x))]
        position_interval = seq(from = -K, to = K, by = 1L)
        frequency_interval = position_interval
        elements = build_f_g_Ff_Fg(type, sigma)
        elements = build_Sf_Sg(elements, lambda, position_interval, frequency_interval)
        ## Create data to plot
        if(("closed" %in% to_plot) & (!is.null(elements$Sf_closed))) {
          Sf_closed = data.frame(x = x, y = elements$Sf_closed(x), method = "closed")
        } else {
          Sf_closed = data.frame(x = x, y = rep(NA, length(x)), method = "closed")
        }
        if("direct" %in% to_plot) {
          Sf_approx_direct = data.frame(x = x, y = elements$Sf_approx_direct(x), method = "direct")
        } else {
          Sf_approx_direct = data.frame(x = x, y = rep(NA, length(x)), method = "direct")
        }
        if("Fourier" %in% to_plot) {
          Sf_approx_Fourier = data.frame(x = x, y = elements$Sf_approx_Fourier(x), method = "Fourier")
        } else {
          Sf_approx_Fourier = data.frame(x = x, y = rep(NA, length(x)), method = "Fourier")
        }
        data = rbind(Sf_approx_direct, Sf_approx_Fourier, Sf_closed)
        idx_remove = which(is.na(data$y))
        if(length(idx_remove) > 0) {
          data = data[-idx_remove,]
        }

        get_title = function(type, sigma, lambda) {
          # lambda as a function of sigma
          get_lambda_as_fn_of_sigma = function(sigma, lambda) {
            if(sigma/lambda > 1) {
              fsigma = paste0("σ/", round(sigma/lambda))
            } else if(sigma/lambda == 1) {
              fsigma = "σ"
            } else {
              fsigma = paste0(round(lambda/sigma), "σ")
            }
            paste0("λ=", fsigma)
          }
          # sigma as a function of sigma0
          get_sigma_as_fn_of_sigma0 = function(type, sigma) {
            sigma0 = sigma_such_as_Fourier_tranform_sums_to_one_func(type)
            if(sigma0/sigma > 1) {
              fsigma = paste0("σ₀/", round(sigma0/sigma))
            } else if(sigma0/sigma == 1) {
              fsigma = "σ₀"
            } else {
              fsigma = paste0(round(sigma/sigma0), "σ₀")
            }
            paste0("σ=", fsigma)
          }
          if(get_sigma_as_fn_of_sigma0(type, sigma) == "σ=σ₀") {
            title = paste0(type, " (", get_lambda_as_fn_of_sigma(sigma, lambda), ")")
          } else {
            title = paste0(type, " (", get_sigma_as_fn_of_sigma0(type, sigma), ", ", get_lambda_as_fn_of_sigma(sigma, lambda), ")")
          }
          return(title)
        }
        my_title = get_title(type, sigma, lambda) # (only for debug purpose in ggplot, not shown in the final output)

        break_func = function(lambda, add_lambda_over_two = FALSE) {
          if(add_lambda_over_two) {
            return(scale_x_continuous(
              breaks = c(-lambda, -lambda/2, 0, lambda/2, lambda),
              labels = c("-λ", "-λ/2", "0", "λ/2", "λ")))
          } else {
            return(scale_x_continuous(
              breaks = c(-lambda, 0, lambda),
              minor_breaks = c(-lambda/2, lambda/2),
              labels = c("-λ", "0", "λ")))
          }
        }

        xlab_name = "x"
        ylab_name = expression("S"[λ]*"f"[σ]*"(x)")

        cust_format = function(x) {
          out = gsub("e-0", "e-", format(x,digits=1,nsmall=2))
          if((abs(as.numeric(out)) < 1e-2) & (!grepl("e", out))) {
            out = gsub("e-0", "e-", format(as.numeric(out), scientific = TRUE))
          }
          if(as.numeric(out) == 0) {
            out = "0"
          }
          if(abs(as.numeric(out)) < 1e-4) {
            out = "0"
          }
          if(grepl("e", out)) {
            out = "0"
          }
          if(as.numeric(out) < 0) { # for negative values
            out = substr(out, 1, nchar(out)-1)
          }
          return(out)
        }

        label_small = cust_format(min(data$y))
        if(label_small != "") {
          if(as.numeric(label_small) < 1e-4 & 1/lambda < 0.1) {
            label_small = ""
          }
        }
        my_breaks = c(min(data$y), 1/lambda, max(data$y))
        my_labels = c(label_small,
                      "1/λ ",
                      cust_format(max(data$y)))
        if(diff(range(data$y)) == 0) {
          my_breaks = c(1/lambda)
          my_labels = c(" 1/λ ")
        }

        p = ggplot(data = data, aes(x=x,y=y,color=method,linetype=method)) +
          geom_line() +
          scale_linetype_manual(values = c("closed" = "solid", "Fourier" = "dotted", "direct" = "dashed")) +
          scale_color_manual(values = c("closed"="black", "Fourier"="blue", "direct"="red")) +
          theme_bw() + break_func(lambda, FALSE) +
          ggtitle(my_title) +
          scale_y_continuous(
            breaks = my_breaks,
            minor_breaks = NULL,
            labels = my_labels) +
          ylab(ylab_name) +
          xlab(xlab_name) +
          coord_cartesian(ylim = c(min(data$y), max(data$y)), clip = "off")

        delta_y = max(data$y) - min(data$y)
        delta_x = max(data$x) - min(data$x)

        y_text_lambda = min(data$y) - delta_y*45/100
        my_vjust = 0
        if(delta_y == 0) {
          my_vjust = 6.45
        }

        if(add_lambda) {
          p = p + geom_text(data=data.frame(label = paste0("λ=", format(lambda,digits=1))), mapping = aes(label=label),
                            inherit.aes = FALSE,
                            x = min(data$x) - delta_x*70/100,
                            y = y_text_lambda,
                            hjust = 0,
                            vjust = my_vjust,
                            size = 4)
        }
        return(p)
      }
      type2to_plot = function(type) {
        if(type %in% c("rectangular", "linear", "exponential", "sinc")) {
          # in those cases, we checked on the plots that "direct" and "Fourier" are close
          return(c("closed"))
        } else if(type %in% c("gaussian", "polynomial", "sinc2")) {
          # sinc2 to update once the close form is known
          return(c("Fourier"))
        } else {
          return(c("closed", "direct", "Fourier"))
        }
      }
      types = c("rectangular", "linear", "exponential", "polynomial", "gaussian", "sinc", "sinc2")
      # range of the video to plot
      max_lambda = 30
      min_lambda = 1/3
      nb_lambdas = 100
      type = types[5]
      for(type in types) {
        to_plot = type2to_plot(type)
        sigma = sigma_such_as_Fourier_tranform_sums_to_one_func(type)
        lambdas = exp(seq(from = log(max_lambda), to = log(min_lambda), length.out = nb_lambdas))
        lambdas[51] = 3 # because we would like the static plot to be exactly this one, not 3.091221
        if(type == "linear") {
          lambdas[100] = lambdas[100] - 4e-4 # because linear is constant for 1/3
        }
        output_folder_plots_current = file.path(output_folder_plots, "Sf", type)
        dir.create(output_folder_plots_current, showWarnings = FALSE)
        N = length(lambdas)
        for(k in 1:N) {
          lambda = lambdas[k]
          name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
          cat(paste0(round(lambda, 2), " "))
          p = plot_Sf(type, sigma, lambda, to_plot) +
            theme(legend.position="none") +
            ggtitle(NULL)
          if(type %in% c("gaussian", "polynomial", "sinc2", "rectangular")) {
            p = p +
              scale_color_manual(values = c("closed"="black", "Fourier"="black", "direct"="black")) +
              scale_linetype_manual(values = c("closed" = "solid", "Fourier" = "solid", "direct" = "solid"))
          }
          ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
                 plot = p+force_panelsizes(rows = unit(2, "cm"),cols = unit(2, "cm")),
                 width = 120, height = 120, units = "px", dpi = 85) # larger dpi = larger font

          # Extract specific frames as static plots, and without lambda notation
          if(k %in% c(1, 51, 100)) {
            p2 = plot_Sf(type, sigma, lambda, to_plot, add_lambda = FALSE) +
              theme(legend.position="none") +
              ggtitle(NULL)
            if(type %in% c("gaussian", "polynomial", "sinc2", "rectangular")) {
              p2 = p2 +
                scale_color_manual(values = c("closed"="black", "Fourier"="black", "direct"="black")) +
                scale_linetype_manual(values = c("closed" = "solid", "Fourier" = "solid", "direct" = "solid"))
            }
            ggsave(filename = paste0(output_folder_plots_current, "_", round(lambda,2), ".png"),
                   plot = p2+force_panelsizes(rows = unit(2, "cm"),cols = unit(2, "cm")),
                   width = 120, height = 120, units = "px", dpi = 85) # larger dpi = larger font
          }
        }
      }

      ## Step 2. Create the video for Sf
      save_video_from_png_folder = function(my_folder,
                                            stop_positions_frame_idx = c(1, 51, 100),
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
      for(type in types) {
        my_folder = file.path(output_folder_plots, "Sf", type)
        save_video_from_png_folder(my_folder)
      }

      expect_true(0 == 0)
    }
  }
})

test_that("Sf_closed is correct for the *linear case*, using different formulas", {
  test_formula_basic = function(formula, from_direct = TRUE, boundary_step = NULL, N_tested = 1000) {
    if(is.null(boundary_step)) {
      N = N_tested
      set.seed(1)
      df = data.frame(x = runif(N, -3, 3),
                      lambda = abs(runif(N, -3, 3)),
                      sigma = abs(runif(N, -3, 3)))
    } else {
      step = boundary_step # 1, 1/2, 1/4, 1/8 tested
      x = seq(from = -2, to = 2, by = step)
      lambda = seq(from = step, to = 2, by = step)
      sigma = seq(from = step, to = 2, by = step)
      df = expand.grid(x = x, lambda = lambda, sigma = sigma)
    }
    N = nrow(df)
    Sf_x_closed = rep(NA, N)
    Sf_x_formula = rep(NA, N)
    type = "linear"
    for(i in 1:N) {
      x = df$x[i]
      lambda = df$lambda[i]
      sigma = df$sigma[i]
      if(from_direct) {
        interval = interval=seq(from=-10000L, to = 10000L, by = 1L) # can still fail for very small values
        Sf_x_closed[i] = Sf_approx_direct_func(f_func(type, sigma), lambda, interval)(x)
      } else {
        Sf_x_closed[i] = (Sf_closed_func(type, sigma, lambda))(x)
      }
      Sf_x_formula[i] = formula(x, sigma, lambda)

      if(N > 100) {
        print(i)
      }
      expect_equal(Sf_x_closed[i], Sf_x_formula[i])
    }
  }

  ## Formula 1 (before using the formula for arithmetic series) ----
  formula_1 = function(x, sigma, lambda) {
    M_minus = floor((sigma + x)/lambda)
    M_zero = floor(x/lambda)
    M_plus = floor((sigma - x)/lambda)
    if((-M_minus) > (-M_zero-1)) {
      I_minus = c()
    } else {
      I_minus = (-M_minus):(-M_zero-1)
    }
    if((-M_zero) > M_plus) {
      I_plus = c()
    } else {
      I_plus = (-M_zero):M_plus
    }
    length_I = length(c(I_minus, I_plus))
    (1/sigma)*length_I - (x/sigma^2)*(length(I_plus)-length(I_minus)) - (lambda/sigma^2)*(-sum(I_minus)+sum(I_plus))
  }

  ## Formula 2 (after using the formula for arithmetic series) ----
  formula_2 = function(x, sigma, lambda) {
    M_minus = floor((sigma + x)/lambda)
    M_zero = floor(x/lambda)
    M_plus = floor((sigma - x)/lambda)
    if((-M_minus) > (-M_zero-1)) {
      I_minus = c()
    } else {
      I_minus = (-M_minus):(-M_zero-1)
    }
    if((-M_zero) > M_plus) {
      I_plus = c()
    } else {
      I_plus = (-M_zero):M_plus
    }
    length_I = length(c(I_minus, I_plus))
    (1/sigma)*length_I +
      -(x/sigma^2)*(length(I_plus)-length(I_minus)) +
      - (lambda/(2*sigma^2))*((M_plus - M_zero)*length(I_plus) + (M_minus + M_zero + 1)*length(I_minus))
  }

  ## Formula 3 (before locating only in [-lambda/2, lambda/2)) ----
  formula_3 = function(x, sigma, lambda) {
    M_minus = floor((sigma + x)/lambda)
    M_zero = floor(x/lambda)
    M_plus = floor((sigma - x)/lambda)
    (1/sigma) * (M_minus+M_plus+1) +
      - lambda/(2*sigma^2)*(M_plus*(M_plus+1) - 2*M_zero*(M_zero+1) + M_minus*(M_minus+1)) +
      - (x/sigma^2)*(M_plus - M_minus + 1 + 2*M_zero)
  }

  ## Formula final (also the one used now in the main function) ----
  formula_final = function(x, sigma, lambda) {
    x = tilde_func(x, lambda)
    sigma_tilde = tilde_func(sigma, lambda)
    my_difference = round(abs(x) - abs(sigma_tilde), 10)
    condition = ((my_difference < 0) | ((my_difference == 0) & (sigma_tilde > 0)) | (x == 0))
    (1/lambda)*(1-(sigma_tilde^2/sigma^2)) +
      (condition)*(abs(sigma_tilde)-abs(x))/sigma^2
  }

  ## Test all the formulas (for the linear case)
  # test on a random grid of length N_tested
  test_formula_basic(formula_1, N_tested = 10)
  test_formula_basic(formula_2, N_tested = 10)
  test_formula_basic(formula_3, N_tested = 10)
  test_formula_basic(formula_final, N_tested = 10)
  # test on a grid with step of boundary_step
  test_formula_basic(formula_final, boundary_step = 1/2^0)
})

test_that("Step by step computations to obtain the closed formula for the *linear case*", {
  ## Particularization for x \in [-lambda/2, lambda/2)
  test_in_the_central_interval = function(lambda = 1) {
    sigmas = seq(from = 1/64, to = 10, by = 1/64)
    x = seq(from = -lambda/2, to = lambda/2, by = 1/64)
    x = x[-length(x)]
    df = expand.grid(sigma = sigmas, x = x)
    df$M_plus = floor(round((df$sigma - df$x), 10)/lambda)
    df$M_zero = floor(round(df$x, 10)/lambda)
    df$M_minus = floor(round((df$sigma + df$x), 10)/lambda)
    df$sigma_tilde = sapply(df$sigma, function(sigma){tilde_func(sigma, lambda)})
    df$M_minus_tilde = floor(round(df$sigma_tilde + df$x, 10)/lambda)
    df$M_plus_tilde = floor(round((df$sigma_tilde - df$x), 10)/lambda)
    df$M_minus_tilde_inside = round((df$sigma_tilde + df$x),10)/lambda
    df$M_zero_inside = round(df$x, 10)/lambda
    df$M_plus_tilde_inside = round((df$sigma_tilde - df$x), 10)/lambda

    # We have x in [-lambda/2,lambda/2). Let \tilde{sigma} in [-lambda/2,lambda/2)
    # such that \tilde{sigma} = sigma - i*lambda (for a certain integer i), from
    # which we have: sigma = \tilde{sigma} + i*lambda
    #
    # We have:
    # M_minus = floor((sigma + x)/lambda) = floor((\tilde{sigma} + x)/lambda) + i
    # M_plus = floor((sigma - x)/lambda) = floor((\tilde{sigma} - x)/lambda) + i
    #
    # M_plus - M_minus = floor((\tilde{sigma} - x)/lambda) - floor((\tilde{sigma} + x)/lambda)
    # We let: M_minus_tilde = floor((\tilde{sigma} + x)/lambda)
    #         M_plus_tilde = floor((\tilde{sigma} - x)/lambda)
    # We check this is valid (we need to round with 10 to prevent numerical issues with floor)
    expect_equal(df$M_plus - df$M_minus, df$M_plus_tilde - df$M_minus_tilde)

    # Then we have:
    # M_minus_tilde inside in [-1, 1) so floor in {-1,0}
    # M_zero        inside in [-1/2,1/2) so floor in {-1,0}
    # M_plus_tilde  inside in (-1, 1) so floor in {-1,0}
    expect_true(range(df$M_minus_tilde_inside)[1] >= -1)
    expect_true(range(df$M_minus_tilde_inside)[2] < 1)
    expect_true(range(df$M_zero_inside)[1] >= -0.5)
    expect_true(range(df$M_zero_inside)[2] < 0.5)
    expect_true(range(df$M_plus_tilde_inside)[1] > -1)
    expect_true(range(df$M_plus_tilde_inside)[2] < 1)
    expect_equal(sort(unique(df$M_minus_tilde)), c(-1, 0))
    expect_equal(sort(unique(df$M_zero)), c(-1, 0))
    expect_equal(sort(unique(df$M_plus_tilde)), c(-1, 0))

    # We have M_zero == 0 iff x >= 0, from which we deduce:
    # if x>=0 : 1 + 2*M_zero = 1
    # if x<0  : 1 + 2*M_zero = -1
    # Globally: 1 + 2*M_zero = (-1)^(1_{x<0})
    expect_equal(which(df$M_zero == 0), which(df$x >= 0))
    expect_equal(1+2*df$M_zero, 1-2*(df$x <0))

    ## Terms in x ----
    # We are looking for:
    # S1 = -(x/sigma^2)*(M_plus - M_minus + 1 + 2*M_zero)
    # We have at first:
    # S1 = -(x/sigma^2)*(M_plus_tilde - M_minus_tilde + (-1)^(x < 0))
    # Otherwise, we have:
    # M_plus_tilde = floor((sigma_tilde - x))/lambda)
    # M_minus_tilde = floor((sigma_tilde + x)/lambda)
    # S1 = -(x/sigma^2)*(M_plus_tilde - M_minus_tilde + (-1)^(x < 0))

    # [*Case 1*]: In the case of $$x=0$$, we have $$S_1=0$$.

    # [*Case 2*]: In the case of $$\tilde{\sigma}=0$$ and $$x \neq 0$$,
    # we have $$\tilde{M}^{-}=-\mathbf{1}_{x<0}$$ and
    #         $$\tilde{M}^{+}=-\mathbf{1}_{x>0}$$,
    # so $$M^{+}-M^{-}+1+2M^0 = -\mathbf{1}_{x>0}-\mathbf{1}_{x<0}+1 = 0$$
    # (since $$x \neq 0$$), and finally $$S_1=0$$.
    # In details:
    # If |x|=|σ`| & σ`=0, we get back to x=0, excluded here
    # If |x|<|σ`| & σ`=0, this case is not possible
    # If |x|>|σ`| & σ`=0 & x<0 ==> x<0, M_plus_tilde=0, M_minus_tilde=-1 ==> S1=0
    # If |x|>|σ`| & σ`=0 & x>0 ==> x>0, M_plus_tilde=-1, M_minus_tilde=0 ==> S1=0

    # [*Case 3.1*] |x|<|σ`| & σ`<0 & x<0 ==> x>σ` ==> 0>σ`-x
    # M_plus_tilde = -1
    # M_minus_tilde = -1
    # S1=-(x/σ²)*(-1)=-|x|/σ²

    # [*Case 3.2*] |x|<|σ`| & σ`<0 & x>0 ==> x<-σ` ==> x+σ`<0
    # M_plus_tilde = -1
    # M_minus_tilde = -1
    # S1=-(x/σ²)*(+1)=-|x|/σ²

    # [*Case 3.3*] |x|<|σ`| & σ`>0 & x<0 ==> σ`+x>0
    # M_plus_tilde = 0
    # M_minus_tilde = 0
    # S1=-(x/σ²)*(-1)=-|x|/σ²

    # [*Case 3.4*] |x|<|σ`| & σ`>0 & x>0 ==> σ`-x>0
    # M_plus_tilde = 0
    # M_minus_tilde = 0
    # S1=-(x/σ²)*(+1)=-|x|/σ²

    # [*Case 4.1*] |x|>|σ`| & σ`<0 & x<0 ==> σ`-x>0
    # M_plus_tilde = 0
    # M_minus_tilde = -1
    # S1=0

    # [*Case 4.2*] |x|>|σ`| & σ`<0 & x>0 ==> x+σ`>0
    # M_plus_tilde = -1
    # M_minus_tilde = 0
    # S1=0

    # [*Case 4.3*] |x|>|σ`| & σ`>0 & x<0 ==> x+σ`<0
    # M_plus_tilde = 0
    # M_minus_tilde = -1
    # S1=0

    # [*Case 4.4*] |x|>|σ`| & σ`>0 & x>0 ==> σ`-x<0
    # M_plus_tilde = -1
    # M_minus_tilde = 0
    # S1=0

    # [*Case 5.1*]
    # If (|x|=|σ`| & σ`>0 & x = σ`):
    # M_plus_tilde = 0
    # M_minus_tilde = floor(2σ`/lambda) \in {-1,0}, but all positive so 0
    # In this case: S1 = -(x/sigma^2)*((-1)^(x < 0)) = -|x|/σ²

    # [*Case 5.2*]
    # If (|x|=|σ`| & σ`>0 & x = -σ`)
    # M_plus_tilde = floor(2σ`/lambda) = 0
    # M_minus_tilde = 0
    # In this case: S1 = -(-x/sigma^2) = -|x|/σ²

    # [*Case 6.1*]
    # If (|x|=|σ`| & σ`<0 & x = -σ`)
    # M_plus_tilde = floor(2σ`/lambda) = -1
    # M_minus_tilde = 0
    # In this case: S1 = -(x/sigma^2)*(-1 + 1) = 0

    # [*Case 6.2*]
    # If (|x|=|σ`| & σ`<0 & x = σ`)
    # M_plus_tilde = 0
    # M_minus_tilde = floor(2σ`/lambda) = -1
    # In this case: S1 = -(x/sigma^2)*(0 + 1 + -1) = 0

    # At the end, we have (for all x with [-lambda/2,lambda/2], and all lambda, sigma):
    # If (σ`>0  & |x|=|σ`|), S1 = -|x|/σ²
    # If (σ`<=0 & |x|=|σ`|), S1 = 0
    # If (|x|<|σ`|)        , S1 = -|x|/σ²
    # If (|x|>|σ`|)        , S1 = 0
    # So at the end, we have value independent of lambda:
    # S1 = -(|x|/σ²)*1_{(|x|<|σ`|)|(|x|=|σ`| & σ`>0)}
    df$groundtruth_S1 = (-df$x/df$sigma^2)*(df$M_plus - df$M_minus + 1 + 2*df$M_zero)
    # condition = (|x| < |σ`|)|(σ`>0 & |x| == |σ`|)
    # condition = (my_difference < 0) | ((my_difference == 0) & (df$sigma_tilde > 0))
    df$prediction_S1 = 0
    df$my_difference = round(abs(df$x) - abs(df$sigma_tilde), 10)
    idx_replace = which((df$my_difference < 0) | ((df$my_difference == 0) & (df$sigma_tilde > 0)))
    if(length(idx_replace) > 0) {
      df$prediction_S1[idx_replace] = (-abs(df$x)/df$sigma^2)[idx_replace]
    }
    expect_equal(df$groundtruth_S1, df$prediction_S1)

    ## Constant terms ----
    # Since M0 is in {-1,0}, M0(M0+1)=0
    expect_equal(df$M_zero*(df$M_zero+1), rep(0, nrow(df)))
    df$i = (df$sigma - df$sigma_tilde)/lambda
    expect_equal(df$M_plus, df$M_plus_tilde + df$i)
    expect_equal(df$M_minus, df$M_minus_tilde + df$i)
    expect_equal(df$M_plus_tilde*(df$M_plus_tilde+1), rep(0, nrow(df)))
    expect_equal(df$M_minus_tilde*(df$M_minus_tilde+1), rep(0, nrow(df)))

    # M_plus*(M_plus+1) = (M_plus_tilde+i)(M_plus_tilde+1+i)
    #                   = M_plus_tilde(M_plus_tilde+1)+2iM_plus_tilde+i+i^2
    #                   = 2iM_plus_tilde+i+i^2
    expect_equal(df$M_plus*(df$M_plus+1), 2*df$M_plus_tilde*df$i + df$i^2+ df$i)
    expect_equal(df$M_minus*(df$M_minus+1), 2*df$M_minus_tilde*df$i + df$i^2+ df$i)

    # Middle term interior:
    expect_equal(df$M_plus*(df$M_plus+1)+
                   -2*df$M_zero*(df$M_zero+1)+
                   df$M_minus*(df$M_minus+1),
                 2*df$M_plus_tilde*df$i + 2*df$i^2+ 2*df$i +
                   2*df$M_minus_tilde*df$i)

    # So the middle term is:
    df$groundtruth_S2 = -lambda/(2*df$sigma^2)*(
      df$M_plus*(df$M_plus+1)-2*df$M_zero*(df$M_zero+1)+df$M_minus*(df$M_minus+1)
    )
    expect_equal(df$groundtruth_S2,
                 -lambda/(2*df$sigma^2)*(2*df$M_plus_tilde*df$i + 2*df$i^2+ 2*df$i +
                                           2*df$M_minus_tilde*df$i))
    expect_equal(df$groundtruth_S2,
                 -df$i*lambda/(df$sigma^2)*(df$M_plus_tilde+df$M_minus_tilde+df$i+1))

    ## For the left term
    df$groundtruth_S3 = (1/df$sigma)*(df$M_plus+df$M_minus+1)
    # M_minus = M_minus_tilde+i
    # M_plus = M_plus_tilde+i
    expect_equal(df$groundtruth_S3,
                 (1/df$sigma)*(df$M_plus_tilde+df$M_minus_tilde+2*df$i+1))

    ## For all the constant terms
    df$constant = df$groundtruth_S2 + df$groundtruth_S3
    expect_equal(df$constant, (1/df$sigma)*(
      (-df$i*lambda/df$sigma)*(1+df$M_plus_tilde+df$M_minus_tilde+df$i) +
        (1+df$M_plus_tilde+df$M_minus_tilde+2*df$i)
    ))
    # we have i = (sigma - sigma_tilde)/lambda
    expect_equal(df$constant, (1/df$sigma)*(
      (-1 + (df$sigma_tilde/df$sigma))*(1+df$M_plus_tilde+df$M_minus_tilde+df$i) +
        (1+df$M_plus_tilde+df$M_minus_tilde+2*df$i)
    ))
    expect_equal(df$constant, (1/df$sigma)*(
      (-1 + (df$sigma_tilde/df$sigma))*(1+df$M_plus_tilde+df$M_minus_tilde+df$i) +
        (1+df$M_plus_tilde+df$M_minus_tilde+df$i) + df$i
    ))
    expect_equal(df$constant, (1/df$sigma)*(
      (df$sigma_tilde/df$sigma)*(1+df$M_plus_tilde+df$M_minus_tilde+df$i) + df$i
    ))
    expect_equal(df$constant, (1/df$sigma)*(
      (df$sigma_tilde/df$sigma)*(1+df$M_plus_tilde+df$M_minus_tilde) + (1+df$sigma_tilde/df$sigma)*df$i
    ))
    expect_equal(df$constant, (1/df$sigma)*(
      (df$sigma_tilde/df$sigma)*(1+df$M_plus_tilde+df$M_minus_tilde) + (1+df$sigma_tilde/df$sigma)*(df$sigma - df$sigma_tilde)/lambda
    ))
    expect_equal(df$constant, (1/df$sigma)*(
      (df$sigma_tilde/df$sigma)*(1+df$M_plus_tilde+df$M_minus_tilde) + (1/lambda)*(1+df$sigma_tilde/df$sigma)*(df$sigma - df$sigma_tilde)
    ))
    expect_equal(df$constant, (1/df$sigma)*(
      (df$sigma_tilde/df$sigma)*(1+df$M_plus_tilde+df$M_minus_tilde) + (1/lambda)*(df$sigma-df$sigma_tilde^2/df$sigma)
    ))
    expect_equal(df$constant,
                 (df$sigma_tilde/df$sigma^2)*(1+df$M_plus_tilde+df$M_minus_tilde) + (1/lambda)*(1-df$sigma_tilde^2/df$sigma^2)
    )

    # Now we compute (1+M_plus_tilde+M_minus_tilde) \in {-1,0,1}
    # We update the condition by adding x==0
    # We obtain: (1+M_plus_tilde+M_minus_tilde) = (2*(df$sigma_tilde >= 0)-1)*(df$condition)
    df$condition = ((df$my_difference < 0) | ((df$my_difference == 0) & (df$sigma_tilde > 0)) | (df$x == 0))
    expect_equal(
      (2*(df$sigma_tilde >= 0)-1)*(df$condition),
      (1+df$M_plus_tilde+df$M_minus_tilde)
    )

    # Detail for M_plus_tilde + M_minus_tilde + 1 ----
    # S_0 = df$M_plus_tilde + df$M_minus_tilde + 1
    df$groundtruth_S0 = df$M_plus_tilde + df$M_minus_tilde + 1
    # M_plus_tilde = floor((sigma_tilde - x))/lambda)
    # M_minus_tilde = floor((sigma_tilde + x)/lambda)
    # [*Case 1* / condition realized]: In the case of $$x=0$$
    # M_plus_tilde = M_minus_tilde = floor(sigma_tilde/lambda) = -1_{sigma_tilde<0}
    # so M_plus_tilde+M_minus_tilde+1= 2*(sigma_tilde >= 0)-1
    expect_true(all(dplyr::pull(dplyr::filter(df, x==0, sigma_tilde>=0), groundtruth_S0) == 1))
    expect_true(all(dplyr::pull(dplyr::filter(df, x==0, sigma_tilde<0), groundtruth_S0) == -1))
    # the rest is done as before (long to detail all the cases)
    # End detail for M_minus_tilde + M_plus_tilde + 1 ----

    expect_equal(df$constant,
                 (df$sigma_tilde/df$sigma^2)*(2*(df$sigma_tilde >= 0)-1)*(df$condition) + (1/lambda)*(1-df$sigma_tilde^2/df$sigma^2)
    )
    expect_equal(df$constant,
                 (df$sigma_tilde/df$sigma^2)*(2*(df$sigma_tilde >= 0)-1)*(df$condition) + (1/lambda)*(1-df$sigma_tilde^2/df$sigma^2)
    )
    # since it's val*(2*(val>=0)-1), it's |val| in all cases, even in 0: (for val=sigma_tilde)
    expect_equal(df$constant,
                 (abs(df$sigma_tilde)/df$sigma^2)*(df$condition) + (1/lambda)*(1-df$sigma_tilde^2/df$sigma^2)
    )
    expect_equal(df$constant, (1/lambda) +
                   (abs(df$sigma_tilde)/df$sigma^2)*(df$condition) - (1/lambda)*(df$sigma_tilde^2/df$sigma^2)
    )
    expect_equal(df$constant, (1/lambda) + abs(df$sigma_tilde)*(
      (1/df$sigma^2)*(df$condition) - abs(df$sigma_tilde)/(lambda*df$sigma^2)
    ))
    expect_equal(df$constant, (1/lambda) + (abs(df$sigma_tilde)/df$sigma^2)*(
      (df$condition) - abs(df$sigma_tilde)/lambda
    ))
    expect_equal(df$constant, (1/lambda) + (abs(df$sigma_tilde)/df$sigma^2)*(
      (df$condition) - abs(df$sigma_tilde)/lambda
    ))

    ## Globally ----
    expect_equal(df$groundtruth_S1, (-abs(df$x)/df$sigma^2)*(df$condition))
    expect_equal(df$constant, (1/lambda) + (abs(df$sigma_tilde)/df$sigma^2)*((df$condition) - abs(df$sigma_tilde)/lambda))

    df$Sf = df$groundtruth_S1+df$constant
    expect_equal(df$Sf, (-abs(df$x)/df$sigma^2)*(df$condition)+(1/lambda) + (abs(df$sigma_tilde)/df$sigma^2)*((df$condition) - abs(df$sigma_tilde)/lambda))
    # -(|x|/σ²) 1_{condition} + (1/λ) + (|σ`|/σ²)(1_{condition} - |σ`|/λ)
    # -(|x|/σ²) 1_{condition} + (1/λ) + (|σ`|/σ²)(1_{condition}) -  (σ`/σ)² (1/λ)
    # 1_{condition} (|σ`| - |x|)/σ²  + (1/λ) (1 - (σ`/σ)²)
    expect_equal(df$Sf,
                 (1/lambda)*(1-(df$sigma_tilde^2/df$sigma^2)) +
                   (df$condition)*(abs(df$sigma_tilde)-abs(df$x))/df$sigma^2
    )
    # 1_{condition} (|σ`| - |x|)/σ²  + (1/λ) (1 - (σ`/σ)²)
    # with: condition = (|x|<|σ`| or (|x|=|σ`| and σ`>0) or (x=0))
  }

  # test_in_the_central_interval(lambda = 4)
  # test_in_the_central_interval(lambda = 2)
  # test_in_the_central_interval(lambda = 1/2^0)
  # test_in_the_central_interval(lambda = 1/2^1)
  test_in_the_central_interval(lambda = 1/2^2)
  test_in_the_central_interval(lambda = 1/2^3)
})

test_that("Sf_closed and Sg_closed are correct for the *rectangular case*", {
  formula_rectangular_f = function(x, sigma, lambda) {
    type = "rectangular"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_rectangular_g = function(x, sigma, lambda) {
    type = "rectangular"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  # below on the grid (not the random samples), points are undefined,
  # in those cases we evaluate the derivative as 0

  type = "rectangular"

  # random tests
  N_tested = 10 # 1000
  test_formula(formula_rectangular_f, N_tested = N_tested,
               type = type, original_function = f_func)
  test_formula(formula_rectangular_g, N_tested = N_tested,
               type = type, original_function = g_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_rectangular_f, boundary_step = boundary_step,
               type = type, original_function = f_func)
  test_formula(formula_rectangular_g, boundary_step = boundary_step,
               type = type, original_function = g_func)
})

test_that("Sf_closed and Sg_closed are correct for the *linear case*", {
  formula_linear_f = function(x, sigma, lambda) {
    type = "linear"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_linear_g = function(x, sigma, lambda) {
    type = "linear"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  # below on the grid (not the random samples), points are undefined,
  # in those cases we evaluate the derivative as 0

  type = "linear"

  # random tests
  N_tested = 10 # 300
  test_formula(formula_linear_f, N_tested = N_tested,
               type = type, original_function = f_func)
  test_formula(formula_linear_g, N_tested = N_tested,
               type = type, original_function = g_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_linear_f, boundary_step = boundary_step,
               type = type, original_function = f_func)
  test_formula(formula_linear_g, boundary_step = boundary_step,
               type = type, original_function = g_func)
})

test_that("Sf_closed and Sg_closed are correct for the *exponential case*", {
  formula_exponential_f = function(x, sigma, lambda) {
    type = "exponential"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_exponential_g = function(x, sigma, lambda) {
    type = "exponential"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  # below on the grid (not the random samples), points are undefined,
  # in those cases we evaluate the derivative as 0

  type = "exponential"

  # random tests
  N_tested = 10 # 1000
  test_formula(formula_exponential_f, N_tested = N_tested,
               type = type, original_function = f_func)
  test_formula(formula_exponential_g, N_tested = N_tested,
               type = type, original_function = g_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_exponential_f, boundary_step = boundary_step,
               type = type, original_function = f_func)
  test_formula(formula_exponential_g, boundary_step = boundary_step,
               type = type, original_function = g_func)
})

test_that("Sf_closed and Sg_closed are correct for the *polynomial case*", {
  formula_polynomial_f = function(x, sigma, lambda) {
    type = "polynomial"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_polynomial_g = function(x, sigma, lambda) {
    type = "polynomial"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  type = "polynomial"

  # random tests
  N_tested = 10 # 1000

  test_formula(formula_polynomial_f, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_polynomial_g, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_polynomial_f, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_polynomial_g, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)
})

test_that("Sf_closed and Sg_closed are correct for the *gaussian case*", {
  formula_gaussian_f = function(x, sigma, lambda) {
    type = "gaussian"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_gaussian_g = function(x, sigma, lambda) {
    type = "gaussian"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  type = "gaussian"

  # random tests
  N_tested = 10 # 1000

  test_formula(formula_gaussian_f, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_gaussian_g, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_gaussian_f, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_gaussian_g, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)
})

test_that("Sf_closed and Sg_closed are correct for the *sinc case*", {
  formula_sinc_f = function(x, sigma, lambda) {
    type = "sinc"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_sinc_g = function(x, sigma, lambda) {
    type = "sinc"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  type = "sinc"

  # random tests
  N_tested = 10 # 1000

  test_formula(formula_sinc_f, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_sinc_g, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_sinc_f, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_sinc_g, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)
})

test_that("Sf_closed and Sg_closed are correct for the *sinc2 case*", {
  formula_sinc2_f = function(x, sigma, lambda) {
    type = "sinc2"
    Sf_closed_func(type, sigma, lambda)(x)
  }
  formula_sinc2_g = function(x, sigma, lambda) {
    type = "sinc2"
    Sg_closed_func(type, sigma, lambda)(x)
  }

  type = "sinc2"

  # random tests
  N_tested = 10 # 1000
  test_formula(formula_sinc2_f, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_sinc2_g, N_tested = N_tested,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)

  # grid test
  boundary_step = 1/2^0 # 1/2^2
  test_formula(formula_sinc2_f, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sf_approx_Fourier_func, Ff_func needs to be plug
               Sf_approx = Sf_approx_Fourier_func)
  test_formula(formula_sinc2_g, boundary_step = boundary_step,
               type = type, original_function = Ff_func, # for Sg_approx_Fourier_func, Ff_func needs to be plug too
               Sf_approx = Sg_approx_Fourier_func)
})
