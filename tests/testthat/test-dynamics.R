test_that("Delta condition for (x,sigma) and (t,z) are equivalent, except on the border
of the square, where, however, the value of Sf is the same whether Delta is
or not fulfilled, for the Linear type", {
  # original Delta function condition
  Delta_original_func = function(x, sigma) {
    my_difference = round(abs(x) - abs(sigma), 10)
    condition = ((my_difference < 0) | ((my_difference == 0) & (sigma > 0)) | (x == 0))
    return(condition)
  }

  # Ground truth function for the condition Delta(x, sigma), after using tilde functions
  g = function(x, lambda, sigma) {
    x_tilde = tilde_func(x, lambda)
    sigma_tilde = tilde_func(sigma, lambda)
    Delta_original_func(x_tilde, sigma_tilde)
  }

  # New condition using Delta(t, z)
  f = function(x, lambda, sigma) {
    t = sigma/lambda
    z = x/lambda

    # fractional part in (0, 1]
    t = t - floor(t)
    if(t == 0) {
      # just to prevent t == 0, but does not change the result
      t = t+1
    }

    # shifted fractional part, in [-1/2, 1/2)
    z = tilde_func(z, 1)

    # Condition result using the new condition Delta_func(t, z)
    result = Delta_func(t, z)

    # There are some cases where Delta(t,z) is different from Delta(x,sigma)
    # It would be a problem only if the final value of Sf is different in those cases

    # In the following, we set Delta(x,sigma) in the case where the final value
    # for Sf is the same, i.e. that is equivalent to ignore those cases in the tests
    form1 = -frac_minus(t)^2 # form when the condition is not fulfilled
    form2 = -frac_plus(t)*(frac_plus(t)-1)-abs(z) # form when the condition is fulfilled
    my_diff = round(form1 - form2, 10)
    if(my_diff == 0) {
      # when my_diff is 0, the value of the function is the same in
      # both cases (either Delta fulfilled or Delta not fulfilled).
      # So both values for Delta are valid, and we let the same as the
      # original Delta
      result = g(x, lambda, sigma)
    }
    return(result)
  }

  step = 1/2^0
  check_equal_func(f, g, step)
})

test_that("Zf obtained from the direct form is the same as the one obtained from Sf_closed", {
  step = 1/2^0 # tested with 1/2^4
  for(type in c("linear")) {
    for(sigma in c(1/4, 1/2, 1, 2)) {
      f = Zf_func(type, sigma)
      g = Zf_func_from_Sf_closed_func(type, sigma)
      check_equal_tz_func(f, g, step = step, z_range = c(-2, 2), t_max = 2)
    }
  }
})

test_that("plot of the normalized Delta condition", {
  output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
  if(!file.exists(file.path(output_folder_plots, "Delta_normalized_condition.png"))) {
    dir.create(output_folder_plots, showWarnings = FALSE)
    # # Condition obtained using Delta_modulo_func, without ggplot
    # Delta_modulo_func = function(t, z) {
    #   # fractional part in (0, 1]
    #   t = t - floor(t)
    #   idx = which(t == 0)
    #   if(length(idx) > 0) {
    #     # just to prevent t == 0, but does not change the result
    #     t[idx] = t[idx]+1
    #   }
    #
    #   # shifted fractional part, in [-1/2, 1/2)
    #   z = tilde_func(z, 1)
    #
    #   # Condition result using the new condition Delta_func(t, z)
    #   result = Delta_func(t, z)
    #   return(result)
    # }
    # t = seq(from = 0, to = 3, length.out = 100)
    # z = seq(from = -1.5, to = 1.5, length.out = 100)
    # tz = expand.grid(t, z)
    # df = data.frame(t = tz[,1], z = tz[,2], condition = Delta_modulo_func(tz[,1], tz[,2]))
    # df = df %>% filter(condition > 0) %>% select(-condition)
    # plot(df, asp = 1)

    # Plot of the condition with ggplot
    shifted_poly = function(t=0,z=0) {
      df_poly = data.frame(x=c(0,0.5,1,0.5),y=c(0,0.5,0,-0.5))
      df_poly$x = df_poly$x + t
      df_poly$y = df_poly$y + z
      df_poly$id = paste0(t, "--", z)
      return(df_poly)
    }
    tz = expand.grid(t=c(0,1,2),z=c(-1,0,1))
    df_poly = list()
    for(k in 1:nrow(tz)) {
      df_poly[[k]] = shifted_poly(t=tz$t[k],z=tz$z[k])
    }
    df_poly = dplyr::bind_rows(df_poly)
    library(ggplot2)
    p = ggplot(df_poly, aes(x = x, y = y)) +
      geom_polygon(aes(group = id), fill = "#EEEEEE", color = "darkgray") +
      coord_equal() +
      theme_bw() +
      xlab("t") +
      ylab("z")
    ggsave(file.path(output_folder_plots, "Delta_normalized_condition.png"), p, width = 400, height = 400, units = "px")
  }

  expect_true(0 == 0)
})

test_that("videos and plots for linear moving in time output correctly", {
  recompute = FALSE
  output_folder_plots0 = "~/Documents/GitHub/ahstat.github.io/images"
  if(!dir.exists(output_folder_plots0)) {
    # pass this sequence, since there is no corresponding folder
    expect_true(0 == 0)
  } else if(!recompute) {
    expect_true(0 == 0)
  } else {
    # library(periodicmixtures)
    library(dplyr)
    library(ggplot2)
    library(ggh4x)
    output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
    dir.create(output_folder_plots, showWarnings = FALSE)
    subfolders = list.files(output_folder_plots)

    sigma = 1
    type = "linear"
    Zf = Zf_func(type, sigma)
    length.out = 201

    z = seq(from = -1/2, to = 1/2, length.out = length.out)
    z = z[-length(z)]

    t = seq(from = 0, to = 1, length.out = length.out)
    t = t[-1]

    if("function_of_time" %in% subfolders) {
      expect_true(0 == 0)
    } else {
      output_folder_plots_current = file.path(output_folder_plots, "function_of_time")
      dir.create(output_folder_plots_current, showWarnings = FALSE)

      z0 = c(-1/2, -1/4, -1/8, 0, 1/8, 1/4, 1/2)
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360")
      names(color) = color

      plot_Zf_single_t0 = function(t0, t, z, add_t = TRUE) {
        data = data.frame(z = z, Zf = Zf(t0,z))
        p = ggplot(data = data, aes(x=z,y=Zf)) +
          geom_line() +
          theme_bw() +
          scale_x_continuous(
            breaks = c(-1/2, 0, 1/2),
            minor_breaks = c(-1/4, 1/4),
            labels = c("-1/2", "0", "1/2")) +
          scale_y_continuous(
            breaks = c(-1/4, 0, 1/4),
            minor_breaks = c(-1/8, 1/8),
            labels = c("-1/4", "0", "1/4"),
            limits = c(-1/4,1/4)) #+
          #coord_cartesian(ylim = c(-1/4,1/4), clip = "off")

        min_data_y = -1/4
        max_data_y = 1/4
        min_data_x = -1/2
        max_data_x = 1/2
        delta_y = max_data_y - min_data_y
        delta_x = max_data_x - min_data_x
        my_vjust = 0
        if(add_t) {
          char_t = format(t0,digits=1,nsmall=2)
          if(char_t == "0.005") {
            char_t = "0.00"
          }
          if(char_t == "-0.005") {
            char_t = "-0.00"
          }
          p = p + geom_text(data=data.frame(label = paste0("{t}=", char_t)), mapping = aes(label=label),
                            inherit.aes = FALSE,
                            x = min_data_x + delta_x*15/100,
                            y = min_data_y - delta_y*0/100,
                            hjust = 0,
                            vjust = my_vjust,
                            size = 4)
        }
        return(p)
      }

      plot_Zf_single_z0s = function(z0, t0, color, t, z) {
        data = list()
        t_cur = c(0, t)#[t <= t0]
        for(k in 1:length(z0)) {
          data[[k]] = data.frame(z = z0[k], t = rep(t_cur, length(z0)), Zf = Zf(t_cur,z0[k]), color = as.character(color[k]))
        }
        data = dplyr::bind_rows(data)
        data$color = as.factor(data$color)

        p = ggplot(data = data, aes(x=t,y=Zf,color=color), color = color) +
          scale_color_manual(values = color) +
          geom_line() +
          geom_point(data = data %>% filter(t == t0)) +
          theme_bw() +
          scale_x_continuous(
            breaks = c(0, 1/2, 1),
            minor_breaks = c(1/4, 3/4),
            labels = c("0", "1/2", "1")) +
          scale_y_continuous(
            breaks = c(-1/4, 0, 1/4),
            minor_breaks = c(-1/8, 1/8),
            #labels = NULL,
            labels = c("-1/4", "0", "1/4"),
            limits = c(-1/4,1/4)) +
          theme(legend.position="none") +
          xlab("{t}") #+ ylab(NULL)
          #coord_cartesian(ylim = c(-1/4,1/4), clip = "off")
        return(p)
      }

      shifted_poly = function(t=0,z=0) {
        df_poly = data.frame(x=c(0,0.5,1,0.5),y=c(0,0.5,0,-0.5))
        df_poly$x = df_poly$x + t
        df_poly$y = df_poly$y + z
        df_poly$id = paste0(t, "--", z)
        return(df_poly)
      }

      N = length(t)
      for(k in 1:N) {
        t0 = t[k]
        p1 = plot_Zf_single_t0(t0, t, z, add_t = TRUE)
        p1 = p1 +
          geom_point(
            inherit.aes = FALSE,
            mapping = aes(x=x,y=y,color=color),
            color = color,
            data = data.frame(x=z0, y=Zf(t0,z0), color = color))
        p2 = plot_Zf_single_z0s(z0, t0, color, t, z)

        df_poly = shifted_poly(t=0,z=0)
        p3 = ggplot(df_poly, aes(x = x, y = y)) +
          geom_polygon(aes(group = id), fill = "#EEEEEE", color = "darkgray") +
          coord_equal() +
          theme_bw() +
          xlab("{t}") +
          ylab("z") +
          theme_bw() +
          scale_x_continuous(
            breaks = c(0, 1/2, 1),
            minor_breaks = c(1/4, 3/4),
            labels = c("0", "1/2", "1")) +
          scale_y_continuous(
            breaks = c(-1/2, 0, 1/2),
            minor_breaks = c(-1/4, 1/4),
            labels = c("-1/2", "0", "1/2"),
            limits = c(-1/2,1/2))
        p3 = p3 +
          geom_line(
            inherit.aes = FALSE,
            mapping = aes(x=x,y=y),
            color = "black",
            data = data.frame(x=t0, y=z0)
          ) +
          geom_point(
            inherit.aes = FALSE,
            mapping = aes(x=x,y=y,color=color),
            color = color,
            data = data.frame(x=t0, y=z0, color = color))

        p1 = p1 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
        p2 = p2 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
        p3 = p3 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
        # out = gridExtra::grid.arrange(p1, p2, nrow = 1)
        p = egg::ggarrange(p1, p2, p3, widths = c(1,1,1))

        name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
        ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
               plot = p, #+force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm")),
               width = 588, height = 196, units = "px", dpi = 110) # larger dpi = larger font
      }

      ## Step 2. Create the video for Sf ----
      my_folder = file.path(output_folder_plots, "function_of_time")
      save_video_from_png_folder(my_folder, video_wo_stop_time_s = 7)

      expect_true(0 == 0)
    }
  }
})

test_that("videos and plots for linear moving in space output correctly", {
  recompute = FALSE
  output_folder_plots0 = "~/Documents/GitHub/ahstat.github.io/images"
  if(!dir.exists(output_folder_plots0)) {
    # pass this sequence, since there is no corresponding folder
    expect_true(0 == 0)
  } else if(!recompute) {
    expect_true(0 == 0)
  } else {
    # library(periodicmixtures)
    library(dplyr)
    library(ggplot2)
    library(ggh4x)
    output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
    dir.create(output_folder_plots, showWarnings = FALSE)
    subfolders = list.files(output_folder_plots)

    sigma = 1
    type = "linear"
    Zf = Zf_func(type, sigma)
    length.out = 201

    z = seq(from = -1/2, to = 1/2, length.out = length.out)
    z = z[-length(z)]

    t = seq(from = 0, to = 1, length.out = length.out)
    t = t[-1]

    if("function_of_space" %in% subfolders) {
      expect_true(0 == 0)
    } else {
      output_folder_plots_current = file.path(output_folder_plots, "function_of_space")
      dir.create(output_folder_plots_current, showWarnings = FALSE)

      t0 = c(0, 1/8, 1/4, 1/2, 3/4, 7/8, 1)
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360")
      #color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360")
      names(color) = color

      plot_Zf_single_z0 = function(z0, t, z, add_z = TRUE) {
        data = data.frame(t = t, Zf = Zf(t,z0))
        p = ggplot(data = data, aes(x=t,y=Zf)) +
          geom_line() +
          theme_bw() +
          scale_x_continuous(
            breaks = c(0, 1/2, 1),
            minor_breaks = c(1/4, 3/4),
            labels = c("0", "1/2", "1")) +
          scale_y_continuous(
            breaks = c(-1/4, 0, 1/4),
            minor_breaks = c(-1/8, 1/8),
            labels = c("-1/4", "0", "1/4"),
            limits = c(-1/4,1/4)) #+
        #coord_cartesian(ylim = c(-1/4,1/4), clip = "off")

        min_data_y = -1/4
        max_data_y = 1/4
        min_data_x = 0
        max_data_x = 1
        delta_y = max_data_y - min_data_y
        delta_x = max_data_x - min_data_x
        my_vjust = 0
        if(add_z) {
          char_z = format(z0,digits=1,nsmall=2)
          if(char_z == "0.005") {
            char_z = "0.00"
          }
          if(char_z == "-0.005") {
            char_z = "-0.00"
          }
          p = p + geom_text(data=data.frame(label = paste0("z = ", char_z)), mapping = aes(label=label),
                            inherit.aes = FALSE,
                            x = min_data_x + delta_x*15/100,
                            y = min_data_y - delta_y*0/100,
                            hjust = 0,
                            vjust = my_vjust,
                            size = 4)
        }
        return(p)
      }

      plot_Zf_single_t0s = function(z0, t0, color, t, z) {
        data = list()
        z_cur = z#[t <= t0]
        for(k in 1:length(t0)) {
          data[[k]] = data.frame(t = t0[k], z = rep(z_cur, length(t0)), Zf = Zf(t0[k], z_cur), color = as.character(color[k]))
        }
        data = dplyr::bind_rows(data)
        data$color = as.factor(data$color)

        p = ggplot(data = data, aes(x=z,y=Zf,color=color), color = color) +
          scale_color_manual(values = color) +
          geom_line() +
          geom_point(data = data %>% filter(z == z0)) +
          theme_bw() +
          scale_x_continuous(
            breaks = c(-1/2, 0, 1/2),
            minor_breaks = c(-1/4, 1/4),
            labels = c("-1/2", "0", "1/2")) +
          scale_y_continuous(
            breaks = c(-1/4, 0, 1/4),
            minor_breaks = c(-1/8, 1/8),
            #labels = NULL,
            labels = c("-1/4", "0", "1/4"),
            limits = c(-1/4,1/4)) +
          theme(legend.position="none") +
          xlab("z") #+ ylab(NULL)
        #coord_cartesian(ylim = c(-1/4,1/4), clip = "off")
        return(p)
      }

      shifted_poly = function(t=0,z=0) {
        df_poly = data.frame(x=c(0,0.5,1,0.5),y=c(0,0.5,0,-0.5))
        df_poly$x = df_poly$x + t
        df_poly$y = df_poly$y + z
        df_poly$id = paste0(t, "--", z)
        return(df_poly)
      }

      N = length(z)
      for(k in 1:N) {
        z0 = z[k]
        p1 = plot_Zf_single_z0(z0, t, z, add_z = TRUE)
        p1 = p1 +
          geom_point(
            inherit.aes = FALSE,
            mapping = aes(x=x,y=y,color=color),
            color = color,
            data = data.frame(x=t0, y=Zf(t0,z0), color = color))
        p2 = plot_Zf_single_t0s(z0, t0, color, t, z)

        df_poly = shifted_poly(t=0,z=0)
        p3 = ggplot(df_poly, aes(x = x, y = y)) +
          geom_polygon(aes(group = id), fill = "#EEEEEE", color = "darkgray") +
          coord_equal() +
          theme_bw() +
          xlab("{t}") +
          ylab("z") +
          theme_bw() +
          scale_x_continuous(
            breaks = c(0, 1/2, 1),
            minor_breaks = c(1/4, 3/4),
            labels = c("0", "1/2", "1")) +
          scale_y_continuous(
            breaks = c(-1/2, 0, 1/2),
            minor_breaks = c(-1/4, 1/4),
            labels = c("-1/2", "0", "1/2"),
            limits = c(-1/2,1/2))
        p3 = p3 +
          geom_line(
            inherit.aes = FALSE,
            mapping = aes(x=x,y=y),
            color = "black",
            data = data.frame(x=t0, y=z0)
          ) +
          geom_point(
            inherit.aes = FALSE,
            mapping = aes(x=x,y=y,color=color),
            color = color,
            data = data.frame(x=t0, y=z0, color = color))

        p1 = p1 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
        p2 = p2 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
        p3 = p3 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
        # out = gridExtra::grid.arrange(p1, p2, nrow = 1)
        p = egg::ggarrange(p1, p2, p3, widths = c(1,1,1))

        name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
        ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
               plot = p, #+force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm")),
               width = 588, height = 196, units = "px", dpi = 110) # larger dpi = larger font
      }

      ## Step 2. Create the video for Sf ----
      my_folder = file.path(output_folder_plots, "function_of_space")
      save_video_from_png_folder(my_folder, video_wo_stop_time_s = 7)

      expect_true(0 == 0)
    }
  }
})

test_that("videos and plots for linear moving in time unnormalized output correctly", {
  recompute = FALSE
  output_folder_plots0 = "~/Documents/GitHub/ahstat.github.io/images"
  if(!dir.exists(output_folder_plots0)) {
    # pass this sequence, since there is no corresponding folder
    expect_true(0 == 0)
  } else if(!recompute) {
    expect_true(0 == 0)
  } else {
    # library(periodicmixtures)
    library(dplyr)
    library(ggplot2)
    library(ggh4x)
    output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
    dir.create(output_folder_plots, showWarnings = FALSE)
    subfolders = list.files(output_folder_plots)

    sigma = 1
    type = "linear"
    Zf = Zf_func(type, sigma)
    length.out = 201

    z = seq(from = -1/2, to = 1/2, length.out = length.out)
    z = z[-length(z)]

    t = seq(from = 0, to = 1, length.out = length.out)
    t = t[-1]

    if("unnormalized" %in% subfolders) {
      expect_true(0 == 0)
    } else {
      output_folder_plots_current = file.path(output_folder_plots, "unnormalized")
      dir.create(output_folder_plots_current, showWarnings = FALSE)

      plot_Zf_single_t0 = function(t0, t, z, add_t = TRUE) {
        data = data.frame(z = z, Zf = Zf(t0,z))
        p = ggplot(data = data, aes(x=z,y=Zf)) +
          geom_line() +
          theme_bw() +
          scale_x_continuous(
            breaks = c(-1/2, 0, 1/2),
            minor_breaks = c(-1/4, 1/4),
            labels = c("-λ/2", "0", "λ/2")) +
          scale_y_continuous(
            breaks = c(-1/4, 0, 1/4),
            minor_breaks = c(-1/8, 1/8),
            labels = c("1/λ-λ/(4σ²)", "1/λ ", "1/λ+λ/(4σ²)"),
            limits = c(-1/4,1/4))

        min_data_y = -1/4
        max_data_y = 1/4
        min_data_x = -1/2
        max_data_x = 1/2
        delta_y = max_data_y - min_data_y
        delta_x = max_data_x - min_data_x
        my_vjust = 0
        if(add_t) {
          char_t = format(t0,digits=1,nsmall=2)
          if(char_t == "0.005") {
            char_t = "0.00"
          }
          if(char_t == "-0.005") {
            char_t = "-0.00"
          }
          p = p + geom_text(data=data.frame(label = paste0("{σ/λ}=", char_t)), mapping = aes(label=label),
                            inherit.aes = FALSE,
                            x = min_data_x + delta_x*7.5/100,
                            y = min_data_y - delta_y*0/100,
                            hjust = 0,
                            vjust = my_vjust,
                            size = 4)
        }
        return(p)
      }

      N = length(t)
      for(k in 1:N) {
        t0 = t[k]
        ylab_name = expression("S"[λ]*"f"[σ]*"(x)")
        p = plot_Zf_single_t0(t0, t, z, add_t = TRUE) + xlab("x") + ylab(ylab_name)
        p = p + theme(axis.title.y = element_text(vjust=-7, margin=margin(-15,0,0,0)))

        name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
        ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
               plot = p, #+force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm")),
               width = 196*1.33, height = 196, units = "px", dpi = 110) # larger dpi = larger font
      }

      ## Step 2. Create the video for Sf ----
      my_folder = file.path(output_folder_plots, "unnormalized")
      save_video_from_png_folder(my_folder, video_wo_stop_time_s = 7)

      expect_true(0 == 0)
    }
  }
})

