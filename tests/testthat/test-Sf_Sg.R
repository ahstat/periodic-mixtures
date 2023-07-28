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
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("linear", eps = 7e-3, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("exponential", eps = 2e-3, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("gaussian", eps = 2e-15, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("polynomial", eps = 2e-3, N=1, length.out=1e2, K = 1000L, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc", eps = 2e-2, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_direct_and_Sf_approx_Fourier_gives_same("sinc", eps = 8e-4, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)

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
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("linear", eps = 2, N=1, length.out=1e2, verbose = FALSE), NA) # large Gibbs effect
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("exponential", eps = 0.6, N=1, length.out=1e2, verbose = FALSE), NA) # large Gibbs effect
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("gaussian", eps = 2e-14, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("polynomial", eps = 2e-9, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_direct_and_Sg_approx_Fourier_gives_same("sinc", eps = 6e-4, K = 1000L, N=1, length.out=1e2, verbose = FALSE), NA)

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
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("linear", eps = 4e-16, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("exponential", eps = 9e-11, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 3e-15, approx_is_Fourier = TRUE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sf_approx_and_Sf_closed_gives_same("sinc", eps = 2e-1, K = 1000L, approx_is_Fourier = FALSE, N=1, length.out=1e2, verbose = FALSE), NA)

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
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("linear", eps = 3e-16, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("exponential", eps = 2e-13, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc", eps = 4e-13, approx_is_Fourier = TRUE, N=1, length.out=1e2, verbose = FALSE), NA)
  expect_error(test_Sg_approx_and_Sg_closed_gives_same("sinc", eps = 8e-4, K = 1000L, approx_is_Fourier = FALSE, N=1, length.out=1e2, verbose = FALSE), NA)

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
          # gsub("e-0", "e-", format(x,digits=1,nsmall=2))
          out = gsub("e-0", "e-", format(x,digits=1,nsmall=2))
          if((as.numeric(out) < 1e-2) & !grepl("e", out)) {
            out = gsub("e-0", "e-", format(as.numeric(out), scientific = TRUE))
          }
          if(as.numeric(out) == 0) {
            out = "0"
          }
          if(as.numeric(out) < 1e-4) {
            out = "0"
          }
          if(grepl("e", out)) {
            out = ""
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
        if(type %in% c("exponential", "linear", "sinc")) {
          # in those cases, we checked on the plots that "direct" and "Fourier" are close
          return(c("closed"))
        } else if(type %in% c("gaussian", "polynomial")) {
          return(c("Fourier"))
        } else {
          return(c("closed", "direct", "Fourier"))
        }
      }
      types = c("linear", "exponential", "polynomial", "gaussian", "sinc")
      # range of the video to plot
      max_lambda = 30
      min_lambda = 1/3
      nb_lambdas = 100
      type = types[1]
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
          if(type %in% c("gaussian", "polynomial")) {
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
            if(type %in% c("gaussian", "polynomial")) {
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
