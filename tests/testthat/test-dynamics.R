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
  for(type in c("rectangular", "linear", "exponential", "polynomial", "polynomial_one_shift", "polynomial_two_shifts", "gaussian")) {
    for(sigma in c(1/4, 1/2, 1, 2)) {
      f = Zf_func(type, sigma)
      g = Zf_func_from_Sf_closed_func(type, sigma)
      # t_max = 4 to include two periods for the Rectangular type,
      # and 4 periods for the Linear type
      if(type != "gaussian") {
        check_equal_tz_func(f, g, step = step, z_range = c(-2, 2), t_max = 4)
      } else {
        # there are numerical issues in Zf_func_from_Sf_closed_func for t>3,
        # but it looks fine for t>3 using Zf_func directly
        check_equal_tz_func(f, g, step = step, z_range = c(-2, 2), t_max = 3, tol = 1e-4)
      }

    }
  }
})

test_that("videos and plots for moving in time output correctly", {
  recompute = TRUE
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
    sigma = 1
    type = "gaussian"

    # Begin list of parameters for each type ----
    params = list()

    params[["rectangular"]] = list(
      tmax = 2,
      length.out_z = 501, # smooth the steps
      length.out_t = 201,
      z0 = c(-1/2, -1/4, -1/8, 0, 1/8, 1/4, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_step(),
      breaks = c(-1, 0, 1),
      minor_breaks = c(-1/2, 1/2),
      labels = c("  -1", "0", "   1"),
      limits = c(-1,1),
      char_t_name = "2{t/2}=",
      percent_x = 7/100,

      breaks_z_x = c(0, 1, 2),
      minor_breaks_z_x = c(1/2, 3/2),
      labels_z_x = c("0", "1", "2"),
      breaks_z_y = c(-1, 0, 1),
      minor_breaks_z_y = c(-1/2, 1/2),
      labels_z_y = c("  -1", "0", "   1"),
      limits_z_y = c(-1,1),
      xlab_z = "2{t/2}"
    )

    params[["linear"]] = list(
      tmax = 1,
      length.out_z = 201,
      length.out_t = 201,
      z0 = c(-1/2, -1/4, -1/8, 0, 1/8, 1/4, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_line(),
      breaks = c(-1/4, 0, 1/4),
      minor_breaks = c(-1/8, 1/8),
      labels = c("-1/4", "0", "1/4"),
      limits = c(-1/4,1/4),
      char_t_name = "{t}=",
      percent_x = 15/100,

      breaks_z_x = c(0, 1/2, 1),
      minor_breaks_z_x = c(1/4, 3/4),
      labels_z_x = c("0", "1/2", "1"),
      breaks_z_y = c(-1/4, 0, 1/4),
      minor_breaks_z_y = c(-1/8, 1/8),
      labels_z_y = c("-1/4", "0", "1/4"),
      limits_z_y = c(-1/4,1/4),
      xlab_z = "{t}"
    )

    B = 1/2-1/sqrt(6)
    C = 1/2 - 1/(2*sqrt(3))
    params[["exponential"]] = list(
      tmax = 1,
      length.out_z = 201,
      length.out_t = 201,
      z0 = c(-1/2, -C, -B, 0, B, C, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_line(),
      breaks = c(-1/6, 0, 1/6, 1/3),
      minor_breaks = c(),
      labels = c("-1/6", "0", "1/6", "1/3"),
      limits = c(-1/6,1/3),
      char_t_name = " t=",
      percent_x = 15/100,

      breaks_z_x = c(0, 1/3, 2/3, 1),
      minor_breaks_z_x = c(),
      labels_z_x = c("0", "1/3", "2/3", "1"),
      breaks_z_y = c(-1/6, 0, 1/6, 1/3),
      minor_breaks_z_y = c(),
      labels_z_y = c("-1/6", "0", "1/6", "1/3"),
      limits_z_y = c(-1/6,1/3),
      xlab_z = "t"
    )
    rm(B,C)

    params[["polynomial"]] = list(
      tmax = 2,
      length.out_z = 201,
      length.out_t = 201,
      z0 = c(-1/2, -1/4, -1/6, 0, 1/6, 1/4, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_line(),
      breaks = c(-1, 0, 1, 2), # Zf
      minor_breaks = c(),
      labels = c(-1, 0, 1, 2),
      limits = c(-999,999),
      char_t_name = " t=",
      percent_x = 15/100,

      breaks_z_x = c(0, 1, 2), # time
      minor_breaks_z_x = c(1/2, 3/2),
      labels_z_x = c(0, 1, 2),
      breaks_z_y = c(-1, 0, 1, 2), # Zf
      minor_breaks_z_y = c(),
      labels_z_y = c(-1, 0, 1, 2),
      limits_z_y = c(-999,999),
      xlab_z = "t"
    )

    params[["polynomial_one_shift"]] = list(
      tmax = 2,
      length.out_z = 201,
      length.out_t = 201,
      z0 = c(-1/2, -1/4, -1/6, 0, 1/6, 1/4, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_line(),
      breaks = c(-1, 0, 1), # Zf
      minor_breaks = c(-1/2, 1/2),
      labels = c(-1, 0, 1),
      limits = c(-999,999),
      char_t_name = " t=",
      percent_x = 15/100,

      breaks_z_x = c(0, 1, 2), # time
      minor_breaks_z_x = c(1/2, 3/2),
      labels_z_x = c(0, 1, 2),
      breaks_z_y = c(-1, 0, 1), # Zf
      minor_breaks_z_y = c(-1/2, 1/2),
      labels_z_y = c(-1, 0, 1),
      limits_z_y = c(-999,999),
      xlab_z = "t"
    )

    params[["polynomial_two_shifts"]] = list(
      tmax = 2,
      length.out_z = 201,
      length.out_t = 201,
      z0 = c(-1/2, -1/4, -1/6, 0, 1/6, 1/4, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_line(),
      breaks = c(-1, 0, 1), # Zf
      minor_breaks = c(-1/2, 1/2),
      labels = c(-1, 0, 1),
      limits = c(-999,999),
      char_t_name = " t=",
      percent_x = 15/100,

      breaks_z_x = c(0, 1, 2), # time
      minor_breaks_z_x = c(1/2, 3/2),
      labels_z_x = c(0, 1, 2),
      breaks_z_y = c(-1, 0, 1), # Zf
      minor_breaks_z_y = c(-1/2, 1/2),
      labels_z_y = c(-1, 0, 1),
      limits_z_y = c(-999,999),
      xlab_z = "t"
    )

    params[["gaussian"]] = list(
      tmax = 2,
      length.out_z = 201,
      length.out_t = 201,
      z0 = c(-1/2, -1/4, -1/6, 0, 1/6, 1/4, 1/2),
      color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
      geom_type = geom_line(),
      breaks = c(-1, 0, 1, 2), # Zf
      minor_breaks = c(),
      labels = c(-1, 0, 1, 2),
      limits = c(-999,999),
      char_t_name = " t=",
      percent_x = 15/100,

      breaks_z_x = c(0, 1, 2), # time
      minor_breaks_z_x = c(1/2, 3/2),
      labels_z_x = c(0, 1, 2),
      breaks_z_y = c(-1, 0, 1, 2), # Zf
      minor_breaks_z_y = c(),
      labels_z_y = c(-1, 0, 1, 2),
      limits_z_y = c(-999,999),
      xlab_z = "t"
    )
    # End list of parameters for each type ----

    for(type in c("rectangular", "linear", "exponential", "polynomial", "polynomial_one_shift", "polynomial_two_shifts", "gaussian")) {
      dir.create(file.path(output_folder_plots, type), showWarnings = FALSE)
      subfolders = list.files(file.path(output_folder_plots, type))
      Zf = Zf_func(type, sigma)

      tmax = params[[type]]$tmax
      length.out_z = params[[type]]$length.out_z
      length.out_t = params[[type]]$length.out_t

      z = seq(from = -1/2, to = 1/2, length.out = length.out_z)
      z = z[-length(z)]

      t = seq(from = 0, to = tmax, length.out = length.out_t)
      t = t[-1]

      if("function_of_time" %in% subfolders) {
        expect_true(0 == 0)
      } else {
        output_folder_plots_current = file.path(output_folder_plots, type, "function_of_time")
        dir.create(output_folder_plots_current, showWarnings = FALSE)

        z0 = params[[type]]$z0
        color = params[[type]]$color
        names(color) = color

        plot_Zf_single_t0 = function(t0, t, z, type, add_t = TRUE) {
          data = data.frame(z = z, Zf = Zf(t0,z))
          p = ggplot(data = data, aes(x=z,y=Zf))
          p = p + params[[type]]$geom_type
          p = p + theme_bw() +
            scale_x_continuous(
              breaks = c(-1/2, 0, 1/2),
              minor_breaks = c(-1/4, 1/4),
              labels = c("-1/2", "0", "1/2"))
          min_data_x = -1/2
          max_data_x = 1/2

          breaks = params[[type]]$breaks
          minor_breaks = params[[type]]$minor_breaks
          labels = params[[type]]$labels
          limits = params[[type]]$limits

          p = p +
            scale_y_continuous(
              breaks = breaks,
              minor_breaks = minor_breaks,
              labels = labels,
              limits = limits)
          min_data_y = limits[1]
          max_data_y = limits[2]

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
            # need to shift since the fractional part cannot be equal to 2 (rect) or 1 (lin)
            if(type == "rectangular") {
              if(char_t == "2.00") {
                char_t = "0.00"
              }
            } else if(type == "linear") {
              if(char_t == "1.00") {
                char_t = "0.00"
              }
            }

            percent_y = 0/100

            char_t_name = params[[type]]$char_t_name
            percent_x = params[[type]]$percent_x

            p = p + geom_text(
              data=data.frame(label = paste0(char_t_name, char_t)),
              mapping = aes(label=label),
              inherit.aes = FALSE,
              x = min_data_x + delta_x*percent_x,
              y = min_data_y - delta_y*percent_y,
              hjust = 0,
              vjust = my_vjust,
              size = 4)
          }

          return(p)
        }

        plot_Zf_single_z0s = function(z0, t0, color, t, z, type) {
          data = list()
          if(type != "gaussian") {
            t_cur = c(0, t)#[t <= t0]
          } else {
            t_cur = t
          }

          for(k in 1:length(z0)) {
            data[[k]] = data.frame(z = z0[k],
                                   t = rep(t_cur, length(z0)),
                                   Zf = Zf(t_cur,z0[k]),
                                   color = as.character(color[k]))
          }
          data = dplyr::bind_rows(data)
          data$color = as.factor(data$color)

          p = ggplot(data = data, aes(x=t,y=Zf,color=color), color = color) +
            scale_color_manual(values = color) +
            geom_line() +
            geom_point(data = data %>% filter(t == t0)) +
            theme_bw()

          breaks_z_x = params[[type]]$breaks_z_x
          minor_breaks_z_x = params[[type]]$minor_breaks_z_x
          labels_z_x = params[[type]]$labels_z_x
          breaks_z_y = params[[type]]$breaks_z_y
          minor_breaks_z_y = params[[type]]$minor_breaks_z_y
          labels_z_y = params[[type]]$labels_z_y
          limits_z_y = params[[type]]$limits_z_y
          xlab_z = params[[type]]$xlab_z

          p = p +
            scale_x_continuous(
              breaks = breaks_z_x,
              minor_breaks = minor_breaks_z_x,
              labels = labels_z_x) +
            scale_y_continuous(
              breaks = breaks_z_y,
              minor_breaks = minor_breaks_z_y,
              labels = labels_z_y,
              limits = limits_z_y) +
            theme(legend.position="none") +
            xlab(xlab_z)

          return(p)
        }

        shifted_poly = function(t=0,z=0) {
          df_poly = data.frame(x=c(0,0.5,1,0.5),y=c(0,0.5,0,-0.5))
          df_poly$x = df_poly$x + t
          df_poly$y = df_poly$y + z
          df_poly$id = paste0(t, "--", z)
          return(df_poly)
        }

        condition_plot_func = function(type) {
          if(type == "rectangular") {
            df_poly = shifted_poly(t=0,z=0)
            df_poly$x = 2*df_poly$x # since the function is 2-periodic in time

            p3 = ggplot(df_poly, aes(x = x, y = y)) +
              geom_polygon(aes(group = id), fill = "#EEEEEE") +
              coord_equal() +
              theme_bw() +
              xlab("2{t/2}") +
              ylab("z") +
              theme_bw() +
              scale_x_continuous(
                breaks = c(0, 1, 2),
                minor_breaks = c(1/2, 3/2),
                labels = c("0", "1", "2")) +
              scale_y_continuous(
                breaks = c(-1/2, 0, 1/2),
                minor_breaks = c(-1/4, 1/4),
                labels = c("-1/2", "0", "1/2"),
                limits = c(-1/2,1/2))
            p3 = p3 +
              geom_segment(data.frame(x=0,y=0,xend=1,yend=1/2), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray") +
              geom_segment(data.frame(x=0,y=0,xend=1,yend=-1/2), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray") +
              geom_segment(data.frame(x=1,y=1/2,xend=2,yend=0), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray", linetype = "dashed") +
              geom_segment(data.frame(x=1,y=-1/2,xend=2,yend=0), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray", linetype = "dashed") +
              geom_point(data.frame(x=c(1,1),y=c(-1/2,1/2)), mapping =aes(x=x, y=y), inherit.aes = FALSE, color = "darkgray", fill = "white", shape = 21) +
              geom_point(data.frame(x=c(0),y=c(0)), mapping =aes(x=x, y=y), inherit.aes = FALSE, color = "darkgray", fill = "darkgray", shape = 21)

            p3 = p3 +
              annotate(geom = "text", label = "-1",
                       x = 0.25*c(1,1),
                       y = 0.35*c(-1,1),
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "0",
                       x = 1,
                       y = 0,
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "1",
                       x = (2-0.25)*c(1,1),
                       y = 0.35*c(-1,1),
                       hjust = "center", vjust = "middle")

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
            return(p3)
          } else if(type == "linear") {
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
            return(p3)
          } else if(type == "exponential") {
            curve_zeros_top = function(t) {
              (1-t*acosh(t*sinh(1/t)))/2
            }

            t = seq(from = 0, to = 1.5, length.out = 300)
            df_poly = data.frame(
              x = t,
              ymin = -curve_zeros_top(t),
              ymax = curve_zeros_top(t)
            )
            df_poly[1,]$ymin = 0
            df_poly[1,]$ymax = 0

            p3 = ggplot(df_poly) +
              geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax), fill = "#EEEEEE", color = "darkgray") +
              theme_bw() +
              xlab("t") +
              ylab("z") +
              theme_bw() +
              scale_x_continuous(
                breaks = c(0, 1/3, 2/3, 1),
                minor_breaks = c(),
                labels = c("0", "1/3", "2/3", "1")) +
              scale_y_continuous(
                breaks = c(-1/2, 0, 1/2),
                #breaks = c(-1/2, 0, 1/2-1/(2*sqrt(3)), 1/2),
                minor_breaks = c(-1/4, 1/4),
                labels = c("-1/2", "0", "1/2"),
                #labels = c("-1/2", "0", expression("z"[0](infinity)), "1/2"),
                limits = c(-1/2,1/2)) +
              coord_fixed(xlim=c(0,1))

            p3 = p3 +
              annotate(geom = "text", label = "\u2212",
                       x = 0.5,
                       y = 0.35*c(-1,1), fontface =2,
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "+",
                       x = 0.5,
                       y = 0, fontface =2,
                       hjust = "center", vjust = "middle")

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
            return(p3)
          } else if(grepl("polynomial", type)) {
            curve_zeros_top = function(t) {
              acos(exp(-2*t))/(2*pi)
            }

            t = seq(from = 0, to = 2.5, length.out = 300)
            df_poly = data.frame(
              x = t,
              ymin = -curve_zeros_top(t),
              ymax = curve_zeros_top(t)
            )
            df_poly[1,]$ymin = 0
            df_poly[1,]$ymax = 0

            p3 = ggplot(df_poly) +
              geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax), fill = "#EEEEEE", color = "darkgray") +
              theme_bw() +
              xlab("t") +
              ylab("z") +
              theme_bw() +
              scale_x_continuous(
                breaks = c(0, 1, 2),
                minor_breaks = c(1/2, 3/2),
                labels = c(0, 1, 2)) +
              scale_y_continuous(
                breaks = c(-1/2, 0, 1/2),
                #breaks = c(-1/2, 0, 1/2-1/(2*sqrt(3)), 1/2),
                minor_breaks = c(-1/4, 1/4),
                labels = c("-1/2", "0", "1/2"),
                #labels = c("-1/2", "0", expression("z"[0](infinity)), "1/2"),
                limits = c(-1/2,1/2)) +
              coord_fixed(xlim=c(0,2))

            p3 = p3 +
              annotate(geom = "text", label = "\u2212",
                       x = 1,
                       y = 0.35*c(-1,1), fontface =2,
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "+",
                       x = 1,
                       y = 0, fontface =2,
                       hjust = "center", vjust = "middle")

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
            return(p3)
          } else {
            return(ggplot())
          }
        }

        N = length(t)
        k=1
        for(k in 1:N) {
          t0 = t[k]
          p1 = plot_Zf_single_t0(t0, t, z, type, add_t = TRUE)
          p1 = p1 +
            geom_point(
              inherit.aes = FALSE,
              mapping = aes(x=x,y=y,color=color),
              color = color,
              data = data.frame(x=z0, y=Zf(t0,z0), color = color))
          p2 = plot_Zf_single_z0s(z0, t0, color, t, z, type)
          p3 = condition_plot_func(type)

          if(type == "polynomial") {
            p1 = p1 + coord_cartesian(ylim = c(-1,2))
            p2 = p2 + coord_cartesian(ylim = c(-1,2))
          } else if(type %in% c("polynomial_one_shift", "polynomial_two_shifts")) {
            p1 = p1 + coord_cartesian(ylim = c(-1,1))
            p2 = p2 + coord_cartesian(ylim = c(-1,1))
          } else if(type == "gaussian") {
            p1 = p1 + coord_cartesian(ylim = c(-1,2))
            p2 = p2 + coord_cartesian(ylim = c(-1,2))
          }

          p1 = p1 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
          p2 = p2 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
          p3 = p3 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
          # out = gridExtra::grid.arrange(p1, p2, nrow = 1)
          p = suppressWarnings(egg::ggarrange(p1, p2, p3, widths = c(1,1,1))) # warning of font for - unicode sign

          name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
          ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
                 plot = p, #+force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm")),
                 width = 588, height = 196, units = "px", dpi = 110) # larger dpi = larger font
        }

        ## Step 2. Create the video for Sf ----
        my_folder = file.path(output_folder_plots_current)
        save_video_from_png_folder(my_folder, video_wo_stop_time_s = 7)

        expect_true(0 == 0)
      }
    }
  }
})

test_that("videos and plots for moving in space output correctly", {
  recompute = TRUE
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
    sigma = 1
    type = "polynomial"
    for(type in c("rectangular", "linear", "exponential", "polynomial", "polynomial_one_shift", "polynomial_two_shifts")) {
      dir.create(file.path(output_folder_plots, type), showWarnings = FALSE)
      subfolders = list.files(file.path(output_folder_plots, type))
      Zf = Zf_func(type, sigma)

      # Begin list of parameters for each type ----
      params = list()

      params[["rectangular"]] = list(
        tmax = 2,
        length.out_z = 201,
        length.out_t = 201,
        t0 = c(0, 1/2, 1, 3/2, 2),
        color = c("#154360", "#1ABC9C", "gray", "#FFC300", "#154360"),
        geom_type = geom_step(),
        breaks_z_x = c(0, 1, 2),
        minor_breaks_z_x = c(1/2, 3/2),
        labels_z_x = c("0", "1", "2"),
        breaks_z_y = c(-1, 0, 1),
        minor_breaks_z_y = c(-1/2, 1/2),
        labels_z_y =c("-1", "0", "1"),
        limits_z_y = c(-1,1),

        breaks_t_x = c(-1/2, 0, 1/2),
        minor_breaks_t_x = c(-1/4, 1/4),
        labels_t_x = c("-1/2", "0", "1/2"),
        breaks_t_y = c(-1, 0, 1),
        minor_breaks_t_y = c(-1/2, 1/2),
        labels_t_y = c("-1", "0", "1"),
        limits_t_y = c(-1,1)
      )

      params[["linear"]] = list(
        tmax = 1,
        length.out_z = 201,
        length.out_t = 201,
        t0 = c(0, 1/8, 1/4, 1/2, 3/4, 7/8, 1),
        color = c("#154360", "#1ABC9C", "#FFC300", "gray", "#FFC300", "#1ABC9C", "#154360"),
        geom_type = geom_line(),
        breaks_z_x = c(0, 1/2, 1),
        minor_breaks_z_x = c(1/4, 3/4),
        labels_z_x = c("0", "1/2", "1"),
        breaks_z_y = c(-1/4, 0, 1/4),
        minor_breaks_z_y = c(-1/8, 1/8),
        labels_z_y = c("-1/4", "0", "1/4"),
        limits_z_y = c(-1/4,1/4),

        breaks_t_x = c(-1/2, 0, 1/2),
        minor_breaks_t_x = c(-1/4, 1/4),
        labels_t_x = c("-1/2", "0", "1/2"),
        breaks_t_y = c(-1/4, 0, 1/4),
        minor_breaks_t_y = c(-1/8, 1/8),
        labels_t_y = c("-1/4", "0", "1/4"),
        limits_t_y = c(-1/4,1/4)

      )

      params[["exponential"]] = list(
        tmax = 1,
        length.out_z = 201,
        length.out_t = 201,
        t0 = c(0, 1/8, 1/4, 1),
        color = c("gray", "#1ABC9C", "#FFC300", "#154360"),
        geom_type = geom_line(),
        breaks_z_x = c(0, 1/2, 1),
        minor_breaks_z_x = c(1/4, 3/4),
        labels_z_x = c("0", "1/2", "1"),
        breaks_z_y = c(-1/6, 0, 1/6, 1/3),
        minor_breaks_z_y = c(),
        labels_z_y = c("-1/6", "0", "1/6", "1/3"),
        limits_z_y = c(-1/6,1/3),

        breaks_t_x = c(-1/2, 0, 1/2),
        minor_breaks_t_x = c(-1/4, 1/4),
        labels_t_x = c("-1/2", "0", "1/2"),
        breaks_t_y = c(-1/6, 0, 1/6, 1/3),
        minor_breaks_t_y = c(),
        labels_t_y = c("-1/6", "0", "1/6", "1/3"),
        limits_t_y = c(-1/6,1/3)
      )

      params[["polynomial"]] = list(
        tmax = 2,
        length.out_z = 201,
        length.out_t = 201,
        t0 = c(0, log(2)/2, log(3)/2, 2),
        #t0 = c(0, 1/2, 1, 2),
        color = c("gray", "#1ABC9C", "#FFC300", "#154360"),
        geom_type = geom_line(),
        breaks_z_x = c(0, 1, 2), # t
        minor_breaks_z_x = c(1/2, 3/2),
        labels_z_x = c(0, 1, 2),
        breaks_z_y = c(-1, 0, 1, 2), # Zf
        minor_breaks_z_y = c(),
        labels_z_y = c(-1, 0, 1, 2),
        limits_z_y = c(-1,2),

        breaks_t_x = c(-1/2, 0, 1/2), # z
        minor_breaks_t_x = c(-1/4, 1/4),
        labels_t_x = c("-1/2", "0", "1/2"),
        breaks_t_y = c(-1, 0, 1, 2), # Zf
        minor_breaks_t_y = c(),
        labels_t_y = c(-1, 0, 1, 2),
        limits_t_y = c(-1, 2)
      )

      params[["polynomial_one_shift"]] = list(
        tmax = 2,
        length.out_z = 201, # not to change
        length.out_t = 201,
        # t0 = c(1e-6, 1/4, 1/2, 2),
        t0 = c(1e-6, (1/2)*log(5/3), log(3)/2, 2),
        color = c("gray", "#1ABC9C", "#FFC300", "#154360"),
        geom_type = geom_line(),
        breaks_z_x = c(0, 1, 2), # t
        minor_breaks_z_x = c(1/2, 3/2),
        labels_z_x = c(0, 1, 2),
        breaks_z_y = c(-1, 0, 1), # Zf
        minor_breaks_z_y = c(-1/2, 1/2),
        labels_z_y = c(-1, 0, 1),
        limits_z_y = c(-1,1.001),

        breaks_t_x = c(-1/2, 0, 1/2), # z
        minor_breaks_t_x = c(-1/4, 1/4),
        labels_t_x = c("-1/2", "0", "1/2"),
        breaks_t_y = c(-1, 0, 1), # Zf
        minor_breaks_t_y = c(-1/2, 1/2),
        labels_t_y = c(-1, 0, 1),
        limits_t_y = c(-1, 1.001)
      )


      params[["polynomial_two_shifts"]] = list(
        tmax = 2,
        length.out_z = 201,
        length.out_t = 201,
        t0 = c(1e-3, 1/2, 1, 2),
        color = c("gray", "#1ABC9C", "#FFC300", "#154360"),
        geom_type = geom_line(),
        breaks_z_x = c(0, 1, 2), # t
        minor_breaks_z_x = c(1/2, 3/2),
        labels_z_x = c(0, 1, 2),
        breaks_z_y = c(-1, 0, 1), # Zf
        minor_breaks_z_y = c(-1/2, 1/2),
        labels_z_y = c(-1, 0, 1),
        limits_z_y = c(-1,1),

        breaks_t_x = c(-1/2, 0, 1/2), # z
        minor_breaks_t_x = c(-1/4, 1/4),
        labels_t_x = c("-1/2", "0", "1/2"),
        breaks_t_y = c(-1, 0, 1), # Zf
        minor_breaks_t_y = c(-1/2, 1/2),
        labels_t_y = c(-1, 0, 1),
        limits_t_y = c(-1, 1)
      )
      # End list of parameters for each type ----

      tmax = params[[type]]$tmax
      length.out_z = params[[type]]$length.out_z
      length.out_t = params[[type]]$length.out_t

      z = seq(from = -1/2, to = 1/2, length.out = length.out_z)
      z = z[-length(z)]

      t = seq(from = 0, to = tmax, length.out = length.out_t)
      t = t[-1]

      if(type == "polynomial") {
        t = unique(sort(c(t,
                          seq(from = 0, to = 0.1, length.out = 51)[-1],
                          seq(from = 0, to = 0.02, length.out = 51)[-1],
                          seq(from = 0, to = 0.001, length.out = 51)[-1])))
      }

      if(type == "polynomial_one_shift") {
        t = unique(sort(c(t, seq(from = 0, to = 0.02, length.out = 51)[-1])))
      }

      if("function_of_space" %in% subfolders) {
        expect_true(0 == 0)
      } else {
        output_folder_plots_current = file.path(output_folder_plots, type, "function_of_space")
        dir.create(output_folder_plots_current, showWarnings = FALSE)

        t0 = params[[type]]$t0
        color = params[[type]]$color
        names(color) = color

        plot_Zf_single_z0 = function(z0, t, z, type, add_z = TRUE) {
          data = data.frame(t = t, Zf = Zf(t,z0))

          p = ggplot(data = data, aes(x=t,y=Zf))

          breaks_z_x = params[[type]]$breaks_z_x
          minor_breaks_z_x = params[[type]]$minor_breaks_z_x
          labels_z_x = params[[type]]$labels_z_x
          breaks_z_y = params[[type]]$breaks_z_y
          minor_breaks_z_y = params[[type]]$minor_breaks_z_y
          labels_z_y = params[[type]]$labels_z_y
          limits_z_y = params[[type]]$limits_z_y

          min_data_y = limits_z_y[1]
          max_data_y = limits_z_y[2]
          min_data_x = min(breaks_z_x)
          max_data_x = max(breaks_z_x)

          p = p +
            params[[type]]$geom_type +
            theme_bw() +
            scale_x_continuous(
              breaks = breaks_z_x,
              minor_breaks = minor_breaks_z_x,
              labels = labels_z_x) +
            scale_y_continuous(
              breaks = breaks_z_y,
              minor_breaks = minor_breaks_z_y,
              labels = labels_z_y,
              limits = limits_z_y)

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
            p = p + geom_text(
              data=data.frame(label = paste0("z = ", char_z)),
              mapping = aes(label=label),
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
            data[[k]] = data.frame(
              t = t0[k],
              z = rep(z_cur, length(t0)),
              Zf = Zf(t0[k], z_cur),
              color = as.character(color[k]))
          }
          data = dplyr::bind_rows(data)
          data$color = as.factor(data$color)

          p = ggplot(data = data, aes(x=z,y=Zf,color=color), color = color) +
            scale_color_manual(values = color) +
            geom_line() +
            geom_point(data = data %>% filter(z == z0)) +
            theme_bw() +
            theme(legend.position="none") +
            xlab("z")

          breaks_t_x = params[[type]]$breaks_t_x
          minor_breaks_t_x = params[[type]]$minor_breaks_t_x
          labels_t_x = params[[type]]$labels_t_x
          breaks_t_y = params[[type]]$breaks_t_y
          minor_breaks_t_y = params[[type]]$minor_breaks_t_y
          labels_t_y = params[[type]]$labels_t_y
          limits_t_y = params[[type]]$limits_t_y

          p = p + scale_x_continuous(
            breaks = breaks_t_x,
            minor_breaks = minor_breaks_t_x,
            labels = labels_t_x) +
            scale_y_continuous(
              breaks = breaks_t_y,
              minor_breaks = minor_breaks_t_y,
              labels = labels_t_y,
              limits = limits_t_y)

          return(p)
        }

        shifted_poly = function(t=0,z=0) {
          df_poly = data.frame(x=c(0,0.5,1,0.5),y=c(0,0.5,0,-0.5))
          df_poly$x = df_poly$x + t
          df_poly$y = df_poly$y + z
          df_poly$id = paste0(t, "--", z)
          return(df_poly)
        }

        condition_plot_func = function(type) {
          if(type == "rectangular") {
            df_poly = shifted_poly(t=0,z=0)
            df_poly$x = 2*df_poly$x # since the function is 2-periodic in time

            p3 = ggplot(df_poly, aes(x = x, y = y)) +
              geom_polygon(aes(group = id), fill = "#EEEEEE") +
              coord_equal() +
              theme_bw() +
              xlab("2{t/2}") +
              ylab("z") +
              theme_bw() +
              scale_x_continuous(
                breaks = c(0, 1, 2),
                minor_breaks = c(1/2, 3/2),
                labels = c("0", "1", "2")) +
              scale_y_continuous(
                breaks = c(-1/2, 0, 1/2),
                minor_breaks = c(-1/4, 1/4),
                labels = c("-1/2", "0", "1/2"),
                limits = c(-1/2,1/2))
            p3 = p3 +
              geom_segment(data.frame(x=0,y=0,xend=1,yend=1/2), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray") +
              geom_segment(data.frame(x=0,y=0,xend=1,yend=-1/2), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray") +
              geom_segment(data.frame(x=1,y=1/2,xend=2,yend=0), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray", linetype = "dashed") +
              geom_segment(data.frame(x=1,y=-1/2,xend=2,yend=0), mapping =aes(x=x, y=y, xend=xend, yend=yend), inherit.aes = FALSE, col = "darkgray", linetype = "dashed") +
              geom_point(data.frame(x=c(1,1),y=c(-1/2,1/2)), mapping =aes(x=x, y=y), inherit.aes = FALSE, color = "darkgray", fill = "white", shape = 21) +
              geom_point(data.frame(x=c(0),y=c(0)), mapping =aes(x=x, y=y), inherit.aes = FALSE, color = "darkgray", fill = "darkgray", shape = 21)

            p3 = p3 +
              annotate(geom = "text", label = "-1",
                       x = 0.25*c(1,1),
                       y = 0.35*c(-1,1),
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "0",
                       x = 1,
                       y = 0,
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "1",
                       x = (2-0.25)*c(1,1),
                       y = 0.35*c(-1,1),
                       hjust = "center", vjust = "middle")

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
            return(p3)
          } else if(type == "linear") {
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
            return(p3)
          } else if(type == "exponential") {

            curve_zeros_top = function(t) {
              (1-t*acosh(t*sinh(1/t)))/2
            }

            t = seq(from = 0, to = 1.5, length.out = 300)
            df_poly = data.frame(
              x = t,
              ymin = -curve_zeros_top(t),
              ymax = curve_zeros_top(t)
            )
            df_poly[1,]$ymin = 0
            df_poly[1,]$ymax = 0

            p3 = ggplot(df_poly) +
              geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax), fill = "#EEEEEE", color = "darkgray") +
              theme_bw() +
              xlab("t") +
              ylab("z") +
              theme_bw() +
              scale_x_continuous(
                breaks = c(0, 1/3, 2/3, 1),
                minor_breaks = c(),
                labels = c("0", "1/3", "2/3", "1")) +
              scale_y_continuous(
                breaks = c(-1/2, 0, 1/2),
                #breaks = c(-1/2, 0, 1/2-1/(2*sqrt(3)), 1/2),
                minor_breaks = c(-1/4, 1/4),
                labels = c("-1/2", "0", "1/2"),
                #labels = c("-1/2", "0", expression("z"[0](infinity)), "1/2"),
                limits = c(-1/2,1/2)) +
              coord_fixed(xlim=c(0,1))

            p3 = p3 +
              annotate(geom = "text", label = "\u2212",
                       x = 0.5,
                       y = 0.35*c(-1,1), fontface =2,
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "+",
                       x = 0.5,
                       y = 0, fontface =2,
                       hjust = "center", vjust = "middle")

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
            return(p3)
          } else if(grepl("polynomial", type)) {
            curve_zeros_top = function(t) {
              acos(exp(-2*t))/(2*pi)
            }

            t = seq(from = 0, to = 2.5, length.out = 300)
            df_poly = data.frame(
              x = t,
              ymin = -curve_zeros_top(t),
              ymax = curve_zeros_top(t)
            )
            df_poly[1,]$ymin = 0
            df_poly[1,]$ymax = 0

            p3 = ggplot(df_poly) +
              geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax), fill = "#EEEEEE", color = "darkgray") +
              theme_bw() +
              xlab("t") +
              ylab("z") +
              theme_bw() +
              scale_x_continuous(
                breaks = c(0, 1, 2),
                minor_breaks = c(1/2, 3/2),
                labels = c(0, 1, 2)) +
              scale_y_continuous(
                breaks = c(-1/2, 0, 1/2),
                #breaks = c(-1/2, 0, 1/2-1/(2*sqrt(3)), 1/2),
                minor_breaks = c(-1/4, 1/4),
                labels = c("-1/2", "0", "1/2"),
                #labels = c("-1/2", "0", expression("z"[0](infinity)), "1/2"),
                limits = c(-1/2,1/2)) +
              coord_fixed(xlim=c(0,2))

            p3 = p3 +
              annotate(geom = "text", label = "\u2212",
                       x = 1,
                       y = 0.35*c(-1,1), fontface =2,
                       hjust = "center", vjust = "middle") +
              annotate(geom = "text", label = "+",
                       x = 1,
                       y = 0, fontface =2,
                       hjust = "center", vjust = "middle")

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
            return(p3)
          }
        }

        N = length(z)
        k=1
        for(k in 1:N) {
          z0 = z[k]
          p1 = plot_Zf_single_z0(z0, t, z, type, add_z = TRUE)
          p1 = p1 +
            geom_point(
              inherit.aes = FALSE,
              mapping = aes(x=x,y=y,color=color),
              color = color,
              data = data.frame(x=t0, y=Zf(t0,z0), color = color))
          p2 = plot_Zf_single_t0s(z0, t0, color, t, z)

          df_poly = shifted_poly(t=0,z=0)
          p3 = condition_plot_func(type)

          p1 = p1 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
          p2 = p2 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))
          p3 = p3 + theme(axis.title.y = element_text(vjust=-2, margin=margin(-15,0,0,0)))

          # out = gridExtra::grid.arrange(p1, p2, nrow = 1)
          p = suppressWarnings(egg::ggarrange(p1, p2, p3, widths = c(1,1,1))) # warning of font for - unicode sign

          name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
          ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
                 plot = p,
                 width = 588, height = 196, units = "px", dpi = 110) # larger dpi = larger font
        }

        ## Step 2. Create the video for Sf ----
        my_folder = file.path(output_folder_plots_current)
        save_video_from_png_folder(my_folder, video_wo_stop_time_s = 7)

        expect_true(0 == 0)
      }
    }
  }
})

test_that("videos and plots for moving in time unnormalized output correctly", {
  recompute = TRUE
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
    sigma = 1
    type = "exponential"

    # Begin list of parameters for each type ----
    params = list()

    params[["rectangular"]] = list(
      tmax = 2,
      length.out_z = 501, # smooth the steps
      length.out_t = 201,
      breaks_x = c(-1/2, 0, 1/2),
      minor_breaks_x = c(-1/4, 1/4),
      labels_x = c("-λ/2", "0", "λ/2"),
      percent_x = 7.5/100,
      breaks_y = c(-1, 0, 1),
      minor_breaks_y = c(-1/2, 1/2),
      labels_y = c("1/λ-1/σ", "1/λ ", "1/λ+1/σ"),
      limits_y = c(-1,1),
      text_t = "{σ/(2λ)}="
    )

    params[["linear"]] = list(
      tmax = 1,
      length.out_z = 201,
      length.out_t = 201,
      breaks_x = c(-1/2, 0, 1/2),
      minor_breaks_x = c(-1/4, 1/4),
      labels_x = c("-λ/2", "0", "λ/2"),
      percent_x = 7.5/100,
      breaks_y = c(-1/4, 0, 1/4),
      minor_breaks_y = c(-1/8, 1/8),
      labels_y = c("1/λ-λ/(4σ²)", "1/λ ", "1/λ+λ/(4σ²)"),
      limits_y = c(-1/4,1/4),
      text_t = "{σ/λ}="
    )

    params[["exponential"]] = list(
      tmax = 1.5,
      length.out_z = 201,
      length.out_t = 301,
      breaks_x = c(-1/2, 0, 1/2),
      minor_breaks_x = c(-1/4, 1/4),
      labels_x = c("-λ/2", "0", "λ/2"),
      percent_x = 19/100,
      breaks_y = c(-1/6, 0, 1/3),
      minor_breaks_y = c(1/6),
      labels_y = c("1/λ-λ/(6σ²)", "1/λ ", "1/λ+λ/(3σ²)"),
      limits_y = c(-1/6,1/3),
      text_t = "σ/λ="
    )

    params[["polynomial"]] = list(
      tmax = 1.5,
      length.out_z = 60001,
      length.out_t = 301,
      breaks_x = c(-1/2, 0, 1/2),
      minor_breaks_x = c(-1/4, 1/4),
      labels_x = c("-λ/2", "0", "λ/2"),
      percent_x = 19/100,
      breaks_y = c(-1, 0, 1),
      minor_breaks_y = c(-1/2, 1/2, 3/2, 2),
      labels_y = c("1/λ-(2/λ)×\ne^(-2σ/λ)", "1/λ ", "1/λ+(2/λ)×\ne^(-2σ/λ)"), # ²
      limits_y = c(-1,2),
      text_t = "σ/λ="
    )

    # End list of parameters for each type ----

    for(type in c("rectangular", "linear", "exponential", "polynomial")) {
      dir.create(file.path(output_folder_plots, type), showWarnings = FALSE)
      subfolders = list.files(file.path(output_folder_plots, type))
      Zf = Zf_func(type, sigma)

      tmax = params[[type]]$tmax
      length.out_z = params[[type]]$length.out_z
      length.out_t = params[[type]]$length.out_t

      z = seq(from = -1/2, to = 1/2, length.out = length.out_z)
      z = z[-length(z)]

      t = seq(from = 0, to = tmax, length.out = length.out_t)
      t = t[-1]

      if("unnormalized" %in% subfolders) {
        expect_true(0 == 0)
      } else {
        output_folder_plots_current = file.path(output_folder_plots, type, "unnormalized")
        dir.create(output_folder_plots_current, showWarnings = FALSE)

        plot_Zf_single_t0 = function(t0, t, z, type, add_t = TRUE) {
          data = data.frame(z = z, Zf = Zf(t0,z))
          p = ggplot(data = data, aes(x=z,y=Zf)) +
            geom_line() +
            theme_bw()

          breaks_x = params[[type]]$breaks_x
          minor_breaks_x = params[[type]]$minor_breaks_x
          labels_x = params[[type]]$labels_x
          percent_x = params[[type]]$percent_x
          breaks_y = params[[type]]$breaks_y
          minor_breaks_y = params[[type]]$minor_breaks_y
          labels_y = params[[type]]$labels_y
          limits_y = params[[type]]$limits_y

          p = p +
            scale_x_continuous(
              breaks = breaks_x,
              minor_breaks = minor_breaks_x,
              labels = labels_x) +
            scale_y_continuous(
              breaks = breaks_y,
              minor_breaks = minor_breaks_y,
              labels = labels_y,
              limits = limits_y)
          min_data_y = limits_y[1]
          max_data_y = limits_y[2]
          min_data_x = min(breaks_x)
          max_data_x = max(breaks_x)

          delta_y = max_data_y - min_data_y
          delta_x = max_data_x - min_data_x
          my_vjust = 0
          if(add_t) {
            text_t = params[[type]]$text_t
            if(type == "rectangular") {
              # change from (0,2] to (0,1],
              # also changing 2{t/2} in (0,2] to {t/2} in (0,1]
              # original form:
              # char_t = format(t0,digits=1,nsmall=2)
              # text_t = "2{σ/(2λ)}="
              # updated form:
              char_t = format(t0/2,digits=1,nsmall=2)
            } else {
              char_t = format(t0,digits=1,nsmall=2)
            }
            if(char_t == "0.005") {
              char_t = "0.00"
            }
            if(char_t == "-0.005") {
              char_t = "-0.00"
            }

            p = p + geom_text(data=data.frame(label = paste0(text_t, char_t)), mapping = aes(label=label),
                              inherit.aes = FALSE,
                              x = min_data_x + delta_x*percent_x,
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
          if(type == "polynomial") {
            ylab_name = expression("S"[λ]*"f"[σ]*"(x)       ")
          }
          p = plot_Zf_single_t0(t0, t, z, type, add_t = TRUE) + xlab("x") + ylab(ylab_name)
          p = p + theme(axis.title.y = element_text(vjust=-7, margin=margin(-15,0,0,0)))

          name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
          ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
                 plot = p, #+force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm")),
                 width = 196*1.33, height = 196, units = "px", dpi = 110) # larger dpi = larger font
        }

        ## Step 2. Create the video for Sf ----
        my_folder = file.path(output_folder_plots_current)
        save_video_from_png_folder(my_folder, video_wo_stop_time_s = 7)

        expect_true(0 == 0)
      }
    }
  }
})

test_that("additional plots for the exponential case roots output correctly", {
  recompute = FALSE

  library(ggplot2)
  # We consider z_0 the roots of Z_t f(z) in the exponential case
  # We derived the following positive roots (as a function of t)
  z0_root = function(t) {
    (1-t*acosh(t*sinh(1/t)))/2
  }

  # We first check that the roots are correct ----
  type = "exponential"
  sigma = 1
  t = seq(from = 0.003, to = 1, length.out = 1001)
  my_root = sapply(t, function(t){
    Zf_func_z = function(z) {
      Zf_func(type, sigma)(t, z)
    }
    uniroot(Zf_func_z, c(0, 1/2), tol = 1e-10)$root
  })
  expect_equal(my_root, z0_root(t))

  # We then plot the root function with the asymptotes ----
  if(recompute) {
    t = seq(from = 0, to = 1, length.out = 501)
    z0 = z0_root(t)
    z0_eq0 = -t*log(t)/2
    z0_eqInf = 1/2-1/(2*sqrt(3)) - 1/(180*sqrt(3)*t^2)
    z0_limInf = 1/2-1/(2*sqrt(3))
    df = data.frame(t=t,z0=z0,z0_eq0=z0_eq0,z0_eqInf=z0_eqInf,z0_limInf=z0_limInf)
    df[1,]$z0 = 0

    p = ggplot(df, aes(x=t)) +
      geom_hline(yintercept = z0_limInf, col = "red", lty = 2) +
      geom_line(aes(y=z0), linewidth = 0.6) +
      geom_line(aes(y=z0_eqInf), col = "#FFC300", lty = "dashed", alpha = 1) +
      geom_line(aes(y=z0_eq0), col = "#1ABC9C", lty = "dashed", alpha = 1) +
      coord_cartesian(ylim = c(0, z0_limInf)) +
      theme_bw() +
      scale_y_continuous(
        breaks = c(0, 0.1, z0_limInf),
        minor_breaks = c(0.2),
        labels = c("0", "0.1", "1/2-1/(2√3)")) +
      scale_x_continuous(
        breaks = c(0, 0.5, 1),
        minor_breaks = c(0.25, 0.75),
        labels = c("0", "0.5", "1")) +
      ylab(expression(z[0])) +
      theme(axis.title.y = element_text(vjust=-12, margin=margin(0,0,0,0)))
    print(p)
  }
})

test_that("additional plots for the exponential case derivative roots output correctly", {
  recompute = FALSE
  # The derivative d(Z_t f(z))/dt for the exponential type is called Dt here
  # (valid for t>0 and z in [-1/2, 1/2))

  # We are interested by Dt(t,z)=0 given a fixed z ----
  Dt = function(t, z) {
    u = 1-2*abs(z)
    cosh(u/t)/sinh(1/t) - (u/t)*sinh(u/t)/sinh(1/t) + (cosh(1/t)/t)*cosh(u/t)/(sinh(1/t)^2)-2*t
  }

  # For instance for z=-0.21, there is a single root close to t=0.2 ----
  z = -0.21
  t = seq(from = 0, to = 1, length.out = 301)[-1]
  plot(t, Dt(t, z), type = "l")
  abline(h=0, col = "red")

  # Since it does not seem that a closed-form expression exist, we rely
  # on root approximation.

  # With classic approximation, there are issues for small t, as shown there ----
  zs = seq(from = -0.5, to = 0.5, length.out = 101)
  root = rep(NA, length(zs))
  for(k in 1:length(zs)) {
    # print(k)
    z = zs[k]
    f = function(t) {
      Dt(t, z)
    }
    # https://stackoverflow.com/questions/70200174
    sol = try(uniroot(f, c(3e-3,1), tol = .Machine$double.eps), silent=TRUE)
    root[k] <- if (inherits(sol, "try-error")) NA else sol$root
  }
  # plot(root, zs, type = "l")

  # To better approximate the root, we rely on mpfr, that takes longer time
  if(recompute) {
    library(Rmpfr)
    zs1 = seq(from = -0.5, to = 0.5, length.out = 101)
    zs2 = seq(from = -0.01, to = 0.01, length.out = 101)
    zs3 = seq(from = 0.23, to = 0.24, length.out = 11) # close to z0_prim
    zs3_m = -zs3
    zs = unique(sort(c(zs1, zs2, zs3, zs3_m)))
    root = rep(NA, length(zs))
    for(k in 1:length(zs)) {
      print(paste0(k, "/", length(zs)))
      z = zs[k]
      z = mpfr(z, precBits = 64)
      f = function(t) {
        t = mpfr(t, precBits = 64)
        Dt(t, z)
      }
      sol = try(unirootR(f, c(1e-8,10), tol = .Machine$double.eps), silent=TRUE)
      root[k] <- if (inherits(sol, "try-error")) NA else asNumeric(sol$root)
    }
    root
    plot(root, zs, type = "l", xlab = "t", ylab = "z")
    z0_prim = 1/2-sqrt(1-2*sqrt(2/15))/2
    abline(h = z0_prim, lty = "dashed", col = "red")
    # in the middle, the function t -> Z(t,z) is increasing
    # at top and bottom areas, the function t -> Z(t,z) is decreasing
  }
  expect_true(0 == 0)
})

test_that("additional plots for the polynomial_two_shifts case derivative roots output correctly", {
  recompute = FALSE
  # The derivative d(Z_t f(z))/dt for the polynomial_two_shifts type is called Dt here
  # (valid for t>0 and z in [-1/2, 1/2))
  Dt = function(t, z) {
    A = cos(2*pi*z)
    (2/(-2*A*exp(2*t)+exp(4*t)+1)^2)*(
      exp(4*t)*(2*A^2+2*A-4*t-1)+
        exp(6*t)*(-2*A^2+2*A*t+1)+
        exp(2*t)*(2*A*(t-1)-1)+1)
  }

  # We are interested by Dt(t,z)=0 given a fixed z ----
  # For instance for z=-0.21, there is a single root close to t=0.4 ----
  z = -0.21
  t = seq(from = 0, to = 1, length.out = 301)[-1]
  plot(t, Dt(t, z), type = "l")
  abline(h=0, col = "red")
  # Since it does not seem that a closed-form expression exist, we rely
  # on root approximation.
  # But there is an issue: there are sometimes two roots given a fixed z,
  # that are for z>0 within (0.25, 0.29), but are not found by the root algorithm.
  zs = seq(from = -0.5, to = 0.5, length.out = 101)
  root = rep(NA, length(zs))
  for(k in 1:length(zs)) {
    z = zs[k]
    f = function(t) {
      Dt(t, z)
    }
    sol = try(uniroot(f, c(3e-3,1), tol = .Machine$double.eps), silent=TRUE)
    root[k] <- if (inherits(sol, "try-error")) NA else sol$root
  }
  plot(root, zs, type = "l")
  # there is an issue since there are sometimes multiple roots
  # Here is the detail:
  #           first root                     second root
  # z=0.25    [0.68,0.7]                          /
  # z=0.26    [0.8, 0.85]                     [7.85, 7.95]
  # z=0.27    [1,   1.05]                     [3.5,  4]
  # z=0.28    [1.4, 1.6]                      [2,    2.2]
  #
  # We could rely in finding the multiple roots, but instead it is easier
  # to view the problem as finding z given t (instead of finding t given z)

  # With classic approximation, finding z given t ----
  ts = seq(from = 0, to = 4, length.out = 401)[-1]
  root = rep(NA, length(ts))
  for(k in 1:length(ts)) {
    # print(k)
    t = ts[k]
    f = function(z) {
      Dt(t, z)
    }
    # https://stackoverflow.com/questions/70200174
    sol = try(uniroot(f, c(3e-3,0.5), tol = .Machine$double.eps), silent=TRUE)
    root[k] <- if (inherits(sol, "try-error")) NA else sol$root
  }
  plot(ts, root, type = "l", ylim = c(-max(root), max(root)))
  lines(ts, -root, type = "l")

  # To better approximate the root, we rely on mpfr, that takes longer time
  if(recompute) {
    library(Rmpfr)
    ts1 = seq(from = 0, to = 2, length.out = 101)[-1]
    ts2 = seq(from = 0, to = 0.02, length.out = 11)[-1]
    ts3 = seq(from = 0.24, to = 0.29, length.out = 11) # close to 0.25
    # very slow for convergence back to 0.25
    #ts4 = seq(from = 2, to = 10, by = 1)
    #ts5 = seq(from = 11, to = 100, by = 20)
    ts = unique(sort(c(ts1, ts2, ts3)))
    root = rep(NA, length(ts))
    for(k in 1:length(ts)) {
      print(paste0(k, "/", length(ts)))
      t = ts[k]
      t = mpfr(t, precBits = 64)
      f = function(z) {
        z = mpfr(z, precBits = 64)
        Dt(t, z)
      }
      sol = try(unirootR(f, c(1e-8,0.5), tol = .Machine$double.eps), silent=TRUE)
      root[k] <- if (inherits(sol, "try-error")) NA else asNumeric(sol$root)
    }
    root
    plot(ts, root, type = "l", xlab = "t", ylab = "z")
    z0_prim = 1/4
    abline(h = z0_prim, lty = "dashed", col = "red")
    # in the middle, the function t -> Z(t,z) is increasing
    # at top and bottom areas, the function t -> Z(t,z) is decreasing
  }
  expect_true(0 == 0)
})

test_that("additional plot for the exponential case roots and derivative roots output correctly", {
  type = "exponential"
  output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
  output_folder_plots_current = file.path(output_folder_plots, type, "roots_and_derivative_roots.png")

  if(!file.exists(output_folder_plots_current)) {
    z0 = 1/2-1/(2*sqrt(3)) # root at t=+Inf
    z0_prim = 1/2-sqrt(1-2*sqrt(2/15))/2 # derivative root at t=+Inf

    # for z0_root, there is the closed-form expression ----
    z0_root = function(t) {
      (1-t*acosh(t*sinh(1/t)))/2
    }
    t = seq(from = 0, to = 1.5, length.out = 300)
    df_z0 = data.frame(x = t, ymin = -z0_root(t), ymax = z0_root(t))
    df_z0[1,]$ymin = 0
    df_z0[1,]$ymax = 0

    # for z0_prim_root, we rely on approximations (see one of the previous test) ----
    library(Rmpfr)
    library(dplyr)
    # The derivative d(Z_t f(z))/dt for the exponential type is called Dt here
    # (valid for t>0 and z in [-1/2, 1/2)), and we are interested by Dt(t,z)=0
    # given a fixed z
    Dt = function(t, z) {
      u = 1-2*abs(z)
      cosh(u/t)/sinh(1/t) - (u/t)*sinh(u/t)/sinh(1/t) + (cosh(1/t)/t)*cosh(u/t)/(sinh(1/t)^2)-2*t
    }

    zs1 = seq(from = -0.5, to = 0.5, length.out = 101)
    zs2 = seq(from = -0.01, to = 0.01, length.out = 101)
    zs3 = seq(from = 0.23, to = 0.24, length.out = 11) # close to z0_prim
    zs3_m = -zs3
    zs = unique(sort(c(zs1, zs2, zs3, zs3_m)))
    zs = zs[zs >= 0] # symmetric
    root = rep(NA, length(zs))
    for(k in 1:length(zs)) {
      print(paste0(k, "/", length(zs)))
      z = zs[k]
      z = mpfr(z, precBits = 64)
      f = function(t) {
        t = mpfr(t, precBits = 64)
        Dt(t, z)
      }
      sol = try(unirootR(f, c(1e-8,10), tol = .Machine$double.eps), silent=TRUE)
      root[k] <- if (inherits(sol, "try-error")) NA else asNumeric(sol$root)
    }
    df_z0_prim = data.frame(x = root, ymin = -zs, ymax = zs)
    df_z0_prim[1,]$x = 0

    # additional point to compute
    z = mpfr(z0, precBits = 64)
    f = function(t) {
      t = mpfr(t, precBits = 64)
      Dt(t, z)
    }
    sol = try(unirootR(f, c(1e-8,10), tol = .Machine$double.eps), silent=TRUE)
    t_prim_of_z0 = asNumeric(sol$root) # time when Z(t,z0) changes from decreasing to increasing (to 0 finally at t=+inf)

    # plot itself ----
    df_z0$type = "zeros"
    df_z0_prim$type = "zeros_Dt"
    df = rbind(df_z0, df_z0_prim)

    p = ggplot(df) +
      geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax),
                  data = df %>% filter(type == "zeros"),
                  inherit.aes = FALSE,
                  fill = "#EEEEEE", color = "darkgray") +
      geom_line(aes(x = x, y = ymax),
                data = df %>% filter(type == "zeros_Dt"),
                inherit.aes = FALSE,
                color = "orange") +
      geom_line(aes(x = x, y = ymin),
                data = df %>% filter(type == "zeros_Dt"),
                inherit.aes = FALSE,
                color = "orange") +
      theme_bw() +
      xlab("t") +
      ylab("z") +
      theme_bw() +
      scale_x_continuous(
        breaks = c(0, 0.5, 1, 1.5),
        minor_breaks = c(0.25, 0.75)) +
      scale_y_continuous(
        breaks = c(0, 0.1, 0.2, 0.3),
        #breaks = c(-1/2, 0, 1/2-1/(2*sqrt(3)), 1/2),
        minor_breaks = NULL,
        #labels = c("-1/2", "0", "1/2"),
        #labels = c("-1/2", "0", expression("z"[0](infinity)), "1/2"),
        limits = c(-1/2,1/2)) +
      coord_cartesian(xlim=c(0,1), ylim=c(-0.05,0.3), clip="on")

    p2 = p +
      annotate(geom = "text", label = "\u2212↘",
               x = 0.5,
               y = 0.28, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "\u2212↗",
               x = 0.5,
               y = 0.218, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "+↗",
               x = 0.5,
               y = 0, fontface =2,
               hjust = "center", vjust = "middle") +
      geom_point(data = data.frame(x=1.039, y=z0), aes(x=x,y=y), color = "darkgray")+
      geom_point(data = data.frame(x=1.039, y=z0_prim), aes(x=x,y=y), color="#FFC300")+
      geom_point(data = data.frame(x=t_prim_of_z0, y=z0), aes(x=x,y=y), color="#154360")

    # Gray point:
    # t = +Inf                                   = +Inf
    # z = z0 = 1/2-1/(2*sqrt(3))                 =  0.2113249
    # Z(t,z) = 0 by construction                 =  0

    # Orange point:
    # t = +Inf                                   = +Inf
    # z = z0_prim = 1/2-sqrt(1-2*sqrt(2/15))/2   =  0.2403352
    # Z(t,z) = Zf_func(type,1)(9999, z0_prim)    = -0.03181504 (= 1/3-sqrt(2/15) using 2*z^2-2*z+1/3)

    # Blue point:
    # t = t_prim_of_z0                           =  0.2294541 (no closed form)
    # z = z0 = 1/2-1/(2*sqrt(3))                 =  0.2113249
    # Z(t,z) = Zf_func(type,1)(t_prim_of_z0, z0) = -0.01603669
    # DZ(t,z)/dt = 0 by construction

    # z0 =  # root at t=+Inf
    # z0_prim =  # derivative root at t=+Inf

    # suppress the warning of points outside the plot
    suppressWarnings(ggsave(filename = file.path(output_folder_plots_current),
                            plot = p2,
                            width = 196*2, height = 196, units = "px", dpi = 110)) # larger dpi = larger font
  }
  expect_true(0 == 0)
})

test_that("additional plot for the polynomial case roots and derivative roots output correctly", {
  type = "polynomial"
  output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
  output_folder_plots_current = file.path(output_folder_plots, type, "roots_and_derivative_roots.png")

  if(!file.exists(output_folder_plots_current)) {
    z0 = 1/4 # root at t=+Inf and derivative root at t=+Inf

    # for z0_root, there is the closed-form expression ----
    z0_root = function(t) {
      acos(exp(-2*t))/(2*pi)
    }

    t = seq(from = 0, to = 2.5, length.out = 3000)
    df_z0 = data.frame(x = t, ymin = -z0_root(t), ymax = z0_root(t))

    # for z0_prim_root, we rely on approximations (see one of the previous test) ----
    library(Rmpfr)
    library(dplyr)
    library(ggplot2)

    # Z = (A*exp(2*t) - 1)/(exp(2*t)+exp(-2*t)-2*A)
    # Dt = (2 e^(2 t) (-1 + 2 A e^(2 t) + e^(4 t) - 2 A^2 e^(4 t)))/(1 - 2 A e^(2 t) + e^(4 t))^2

    Dt = function(t, z) { # polynomial
      A = cos(2*pi*z)
      2*exp(2*t)*(-2*A^2*exp(4*t)+2*A*exp(2*t)+exp(4*t)-1)/(-2*A*exp(2*t)+exp(4*t)+1)^2
    }
    # Dt=0
    # 2*exp(2*t)*(-2*A^2*exp(4*t)+2*A*exp(2*t)+exp(4*t)-1)/(-2*A*exp(2*t)+exp(4*t)+1)^2 = 0
    # -2*A^2*exp(4*t)+2*A*exp(2*t)+exp(4*t)-1 = 0
    # -A^2*exp(4*t)+A*exp(2*t)+exp(4*t)/2-1/2 = 0
    # A^2*exp(4*t)-A*exp(2*t)-exp(4*t)/2+1/2 = 0
    # A^2-A*exp(-2*t)-1/2+exp(-4*t)/2 = 0

    root_deriv = function(z) {
      A = cos(2*pi*z)
      t1 = (1/2)*log((A-sqrt(1-A^2))/(2*A^2-1))
      t2 = (1/2)*log((A+sqrt(1-A^2))/(2*A^2-1))
      return(c(t1,t2))
    }

    f = function(z) {
      A = cos(2*pi*z)
      (1/2)*log((A-sqrt(1-A^2))/(2*A^2-1))
    }

    z = seq(from = -0.5, to = 0.5, length.out = 30001)
    res = suppressWarnings(t(sapply(z, root_deriv)))
    res = rbind(data.frame(t = res[,1], z = z, category = 1),
                data.frame(t = res[,2], z = z, category = 2))

    res1 = rbind(res %>% filter(z >= 0, category == 1)) %>% rename(z1 = z) %>% select(-category) %>% filter(!is.na(t))
    res2 = rbind(res %>% filter(z >= 0, category == 2)) %>% rename(z2 = z) %>% select(-category) %>% filter(!is.na(t))
    res3 = rbind(res %>% filter(z < 0, category == 2)) %>% rename(z3 = z) %>% select(-category) %>% filter(!is.na(t))
    res4 = rbind(res %>% filter(z < 0, category == 1)) %>% rename(z4 = z) %>% select(-category) %>% filter(!is.na(t))

    res = merge(res1, res2, by = "t")
    res = merge(res, res3, by = "t")
    res = merge(res, res4, by = "t")

    df_z0_prim0 = data.frame(t = 0, z1 = 0.25, z2 = 0, z3 = 0, z4 = -0.25)
    df_z0_prim = res
    df_z0_prim = rbind(df_z0_prim0, df_z0_prim)

    # plot itself ----
    p = ggplot() +
      geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax),
                  data = df_z0,
                  inherit.aes = FALSE,
                  fill = "#EEEEEE", color = "darkgray") +
      geom_line(aes(x = t, y = z1),
                data = df_z0_prim,
                inherit.aes = FALSE,
                color = "orange") +
      geom_line(aes(x = t, y = z2),
                data = df_z0_prim,
                inherit.aes = FALSE,
                color = "orange") +
      geom_line(aes(x = t, y = z3),
                data = df_z0_prim,
                inherit.aes = FALSE,
                color = "orange") +
      geom_line(aes(x = t, y = z4),
                data = df_z0_prim,
                inherit.aes = FALSE,
                color = "orange") +
      theme_bw() +
      xlab("t") +
      ylab("z") +
      theme_bw() +
      scale_x_continuous(
        breaks = c(0, 1, 2, 4, 6, 8, 10),
        minor_breaks = c(1/2, 3/2)) +
      scale_y_continuous(
        breaks = c(0, 0.1, 0.2, 0.3, 0.4),
        minor_breaks = c(-0.05, 0.05, 0.15, 0.25, 0.35, 0.45),
        limits = c(-1/2,1/2)) +
      coord_cartesian(xlim=c(0,2), ylim=c(-0.05,0.45), clip="on")

    p2 = p +
      annotate(geom = "text", label = "\u2212↘",
               x = 1,
               y = 0.415, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "\u2212↗",
               x = 1,
               y = 0.29, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "+↗",
               x = 1,
               y = 0.165, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "+↘",
               x = 1,
               y = 0, fontface =2,
               hjust = "center", vjust = "middle") +
      geom_point(data = data.frame(x=2.085, y=1/4), aes(x=x,y=y), color = "darkgray")+
      geom_point(data = data.frame(x=2.085, y=3/8), aes(x=x,y=y), color="#FFC300")+
      geom_point(data = data.frame(x=2.085, y=1/8), aes(x=x,y=y), color="#FFC300")+
      geom_point(data = data.frame(x=0, y=1/4), aes(x=x,y=y), color="#FFC300")+
      geom_point(data = data.frame(x=0, y=0), aes(x=x,y=y), color="#FFC300")

    ggsave(filename = file.path(output_folder_plots_current),
           plot = p2,
           width = 196*2, height = 196, units = "px", dpi = 110) # larger dpi = larger font
  }
  expect_true(0 == 0)
})

test_that("additional plot for the polynomial_one_shift case roots and derivative roots output correctly", {
  type = "polynomial_one_shift"
  output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
  output_folder_plots_current = file.path(output_folder_plots, type, "roots_and_derivative_roots.png")

  if(!file.exists(output_folder_plots_current)) {
    z0 = 1/4 # root at t=+Inf and derivative root at t=+Inf

    # for z0_root, there is the closed-form expression ----
    z0_root = function(t) {
      acos(exp(-2*t))/(2*pi)
    }

    t = seq(from = 0, to = 2.5, length.out = 3000)
    df_z0 = data.frame(x = t, ymin = -z0_root(t), ymax = z0_root(t))

    # for z0_prim_root, we rely on approximations (see one of the previous test) ----
    library(Rmpfr)
    library(dplyr)
    library(ggplot2)

    # Z = (1 - 1/exp(2*t))*(A*exp(2*t) - 1)/(exp(2*t)+exp(-2*t)-2*A)
    # Dt = -(2 (-1 + A) e^(2 t) (-1 - 2 e^(2 t) + (1 + 2 A) e^(4 t)))/(1 - 2 A e^(2 t) + e^(4 t))^2

    Dt = function(t, z) { # polynomial_one_shift
      A = cos(2*pi*z)
      -2*(A-1)*exp(2*t)*((2*A+1)*exp(4*t)-2*exp(2*t)-1)/(-2*A*exp(2*t)+exp(4*t)+1)^2
    }
    # Dt=0
    # -2*(A-1)*exp(2*t)*((2*A+1)*exp(4*t)-2*exp(2*t)-1)/(-2*A*exp(2*t)+exp(4*t)+1)^2 = 0
    # (2*A+1)*exp(4*t)-2*exp(2*t)-1 = 0 or A=1
    # A = exp(-4*t)/2+exp(-2*t)-1/2
    # cos(2*pi*z) = exp(-4*t)/2+exp(-2*t)-1/2
    # 2*pi*z = acos(exp(-4*t)/2+exp(-2*t)-1/2)
    # z = acos(exp(-4*t)/2+exp(-2*t)-1/2)/(2*pi)

    root_deriv = function(t) {
      acos(exp(-4*t)/2+exp(-2*t)-1/2)/(2*pi)
    }

    t1 = seq(from = 0, to = 2.5, length.out = 3001)[-1]
    t2 = seq(from = 0, to = 0.05, length.out = 3001)[-1]
    t = sort(unique(c(t1,t2)))
    res = sapply(t, root_deriv)

    df_z0_prim0 = data.frame(t = 0, z1 = 0, z2 = 0)
    df_z0_prim = data.frame(t = t, z1 = -res, z2 = res)
    df_z0_prim = rbind(df_z0_prim0, df_z0_prim)

    # plot itself ----
    p = ggplot() +
      geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax),
                  data = df_z0,
                  inherit.aes = FALSE,
                  fill = "#EEEEEE", color = "darkgray") +
      geom_line(aes(x = t, y = z1),
                data = df_z0_prim,
                inherit.aes = FALSE,
                color = "orange") +
      geom_line(aes(x = t, y = z2),
                data = df_z0_prim,
                inherit.aes = FALSE,
                color = "orange") +
      theme_bw() +
      xlab("t") +
      ylab("z") +
      theme_bw() +
      scale_x_continuous(
        breaks = c(0, 1, 2, 4, 6, 8, 10),
        minor_breaks = c(1/2, 3/2)) +
      scale_y_continuous(
        breaks = c(0, 0.1, 0.2, 0.3, 0.4),
        minor_breaks = c(-0.05, 0.05, 0.15, 0.25, 0.35, 0.45),
        limits = c(-1/2,1/2)) +
      coord_cartesian(xlim=c(0,2), ylim=c(-0.05,0.35), clip="on")


    p2 = p +
      annotate(geom = "text", label = "\u2212↘",
               x = 1,
               y = 0.345, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "\u2212↗",
               x = 1,
               y = 0.265, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "+↗",
               x = 1,
               y = 0, fontface =2,
               hjust = "center", vjust = "middle") +
      geom_point(data = data.frame(x=2.085, y=1/4), aes(x=x,y=y), color = "darkgray")+
      geom_point(data = data.frame(x=2.085, y=1/3), aes(x=x,y=y), color="#FFC300")+
      geom_point(data = data.frame(x=log(1+sqrt(2))/2, y=1/4), aes(x=x,y=y), color="#154360")

    ggsave(filename = file.path(output_folder_plots_current),
           plot = p2,
           width = 196*2, height = 196, units = "px", dpi = 110) # larger dpi = larger font
  }
  expect_true(0 == 0)
})

test_that("additional plot for the polynomial_two_shifts case roots and derivative roots output correctly", {
  type = "polynomial_two_shifts"
  output_folder_plots = "~/Documents/GitHub/ahstat.github.io/images/2023-6-11-Periodic-mixtures/plot3"
  output_folder_plots_current = file.path(output_folder_plots, type, "roots_and_derivative_roots.png")

  if(!file.exists(output_folder_plots_current)) {
    z0 = 1/4 # root at t=+Inf and derivative root at t=+Inf

    # for z0_root, there is the closed-form expression ----
    z0_root = function(t) {
      acos(exp(-2*t))/(2*pi)
    }

    t = seq(from = 0, to = 2.5, length.out = 3000)
    df_z0 = data.frame(x = t, ymin = -z0_root(t), ymax = z0_root(t))

    # for z0_prim_root, we rely on approximations (see one of the previous test) ----
    library(Rmpfr)
    library(dplyr)
    library(ggplot2)

    Dt = function(t, z) { # polynomial_two_shifts
      A = cos(2*pi*z)
      (2/(-2*A*exp(2*t)+exp(4*t)+1)^2)*(
        exp(4*t)*(2*A^2+2*A-4*t-1)+
          exp(6*t)*(-2*A^2+2*A*t+1)+
          exp(2*t)*(2*A*(t-1)-1)+1)
    }

    ts1 = seq(from = 0, to = 2.5, length.out = 101)[-1]
    ts2 = seq(from = 0, to = 0.02, length.out = 11)[-1]
    ts3 = seq(from = 0.02, to = 0.2, length.out = 21)[-1]
    ts4 = seq(from = 0.24, to = 0.29, length.out = 11) # close to 0.25
    # very slow for convergence back to 0.25
    #ts4 = seq(from = 2, to = 10, by = 1)
    #ts5 = seq(from = 11, to = 100, by = 20)
    ts = unique(sort(c(ts1, ts2, ts3, ts4)))
    root = rep(NA, length(ts))
    for(k in 1:length(ts)) {
      print(paste0(k, "/", length(ts)))
      t = ts[k]
      t = mpfr(t, precBits = 64)
      f = function(z) {
        z = mpfr(z, precBits = 64)
        Dt(t, z)
      }
      sol = try(unirootR(f, c(1e-8,0.5), tol = .Machine$double.eps), silent=TRUE)
      root[k] <- if (inherits(sol, "try-error")) NA else asNumeric(sol$root)
    }

    df_z0_prim0 = data.frame(x = 0, ymin = 0, ymax = 0)
    df_z0_prim = data.frame(x = ts, ymin = -root, ymax = root)
    df_z0_prim = rbind(df_z0_prim0, df_z0_prim)

    # additional point to compute: for z0 (quick and simple)
    z = mpfr(z0, precBits = 64)
    f = function(t) {
      t = mpfr(t, precBits = 64)
      Dt(t, z)
    }
    sol = try(unirootR(f, c(1e-8,10), tol = .Machine$double.eps), silent=TRUE)
    t_prim_of_z0 = asNumeric(sol$root) # time when Z(t,z0) changes from decreasing to increasing (to 0 finally at t=+inf)

    # additional point to compute: for z0_prim (difficult)
    t_interval = c(1,3)
    z_interval = c(1e-8,0.5)
    F_tz = Dt
    l = find_maximum_implicit(F_tz, t_interval, z_interval)

    tz0 = maximum_number_to_str(l)
    z0_prim = asNumeric(tz0[2]) # maximum (over z) of the derivative root
    t_prim_of_z0_prim = asNumeric(tz0[1]) # corresponding t

    # We obtain more precisely (verified)
    # z = 0.2808306321626223388
    # t = 1.750595347

    # plot itself ----
    df_z0$type = "zeros"
    df_z0_prim$type = "zeros_Dt"
    df = rbind(df_z0, df_z0_prim)

    p = ggplot(df) +
      geom_ribbon(aes(x = x, ymin= ymin, ymax = ymax),
                  data = df %>% filter(type == "zeros"),
                  inherit.aes = FALSE,
                  fill = "#EEEEEE", color = "darkgray") +
      geom_line(aes(x = x, y = ymax),
                data = df %>% filter(type == "zeros_Dt"),
                inherit.aes = FALSE,
                color = "orange") +
      geom_line(aes(x = x, y = ymin),
                data = df %>% filter(type == "zeros_Dt"),
                inherit.aes = FALSE,
                color = "orange") +
      theme_bw() +
      xlab("t") +
      ylab("z") +
      theme_bw() +
      scale_x_continuous(
        breaks = c(0, 1, 2, 4, 6, 8, 10),
        minor_breaks = c(1/2, 3/2)) +
      scale_y_continuous(
        breaks = c(0, 0.1, 0.2, 0.3, 0.4),
        minor_breaks = c(-0.05, 0.05, 0.15, 0.25, 0.35, 0.45),
        limits = c(-1/2,1/2)) +
      coord_cartesian(xlim=c(0,2), ylim=c(-0.05,0.35), clip="on")

    p2 = p +
      annotate(geom = "text", label = "\u2212↘",
               x = 1,
               y = 0.33, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "\u2212↗",
               x = 1,
               y = 0.25, fontface =2,
               hjust = "center", vjust = "middle") +
      annotate(geom = "text", label = "+↗",
               x = 1,
               y = 0, fontface =2,
               hjust = "center", vjust = "middle") +
      geom_point(data = data.frame(x=2.085, y=z0), aes(x=x,y=y), color = "darkgray")+
      geom_point(data = data.frame(x=t_prim_of_z0_prim, y=z0_prim), aes(x=x,y=y), color="#FFC300")+
      geom_point(data = data.frame(x=t_prim_of_z0, y=z0), aes(x=x,y=y), color="#154360")

    ggsave(filename = file.path(output_folder_plots_current),
           plot = p2,
           width = 196*2, height = 196, units = "px", dpi = 110) # larger dpi = larger font

    # Gray point:
    # t = +Inf                                   = +Inf
    # z = z0 = 1/4                               =  0.25
    # Z(t,z) = 0 by construction                 =  0

    # Orange point:
    # t = t_prim_of_z0_prim                      =  1.750595
    # z = z0_prim                                =  0.2808306
    # Z(t,z) = Zf_func(type,1)(99, z0_prim)      = -0.1925053
    # DZ(t,z)/dt = 0 by construction

    # Blue point:
    # t = t_prim_of_z0                           =  0.6835852
    # z = z0 = 1/4                               =  0.25
    # Z(t,z) = Zf_func(type,1)(t_prim_of_z0, z0) = -0.0949451
    # DZ(t,z)/dt = 0 by construction

    # z0 =  # root at t=+Inf
    # z0_prim =  # derivative root at t=+Inf
  }
  expect_true(0 == 0)
})






debug = FALSE
if(debug) {
  test_that("debug moving in time", {
    library(periodicmixtures)
    library(dplyr)
    library(ggplot2)
    library(ggh4x)
    sigma = 1
    type = "gaussian"

    # Begin list of parameters for each type ----
    params = list()

    params[[type]] = list(
      tmax = 5,
      length.out_z = 201,
      length.out_t = 201,
      geom_type = geom_line()
    )

    output_folder_plots = "~/Documents/debug"

    # End list of parameters for each type ----
    dir.create(file.path(output_folder_plots, type), showWarnings = FALSE, recursive = TRUE)
    subfolders = list.files(file.path(output_folder_plots, type))
    #Zf = Zf_func(type, sigma)
    #Zf = Zf_func_from_Sf_closed_func(type, sigma)

    scaling_term_func = function(t, sigma, type) {
      if(type == "gaussian") {
        return((exp(pi*t^2)/(2*t)))
      }
      if(type == "rectangular") {
        return(sigma)
      } else if(type == "linear") {
        return(sigma) # return(t*sigma)
      } else if(type == "exponential") {
        return(t*sigma)
      } else if(type == "polynomial") {
        return((exp(2*t)/(2*t))*sigma)
      } else if(type == "polynomial_one_shift") {
        return((exp(2*t)-1)/(2*t)*sigma)
      } else if(type == "polynomial_two_shifts") {
        return((exp(2*t)-1-2*t)/(2*t)*sigma)
      }
    }

    Zf =   function(t, z) {
      # normalization in time
      type0 = ifelse(grepl("polynomial", type), "polynomial", type)
      Zf1 = sapply(t, normalization_in_space_and_time_from_Sf_closed_func, z, sigma, type0)
      # normalization in space
      lambda = sigma/t
      normalization = scaling_term_func(t, sigma, type)
      Zf2 = (Zf1 - 1/lambda)*normalization
      Zf2 = as.numeric(Zf2)
      return(Zf2)
    }


    tmax = params[[type]]$tmax
    length.out_z = params[[type]]$length.out_z
    length.out_t = params[[type]]$length.out_t

    z = seq(from = -1/2, to = 1/2, length.out = length.out_z)
    z = z[-length(z)]

    t = seq(from = 0, to = tmax, length.out = length.out_t)
    t = t[-1]

    output_folder_plots_current = file.path(output_folder_plots, type, "function_of_time")
    dir.create(output_folder_plots_current, showWarnings = FALSE)

    plot_Zf_single_t0 = function(t0, t, z, type, add_t = TRUE) {
      data = data.frame(z = z, Zf = Zf(t0,z))
      p = ggplot(data = data, aes(x=z,y=Zf))
      p = p + params[[type]]$geom_type
      p = p + theme_bw()

      return(p)
    }

    N = length(t)
    k=1
    for(k in 1:N) {
      t0 = t[k]
      p = plot_Zf_single_t0(t0, t, z, type, add_t = TRUE)

      p = p + coord_cartesian(ylim = c(-1,1))

      name = paste0(paste(rep("0", nchar(N) - nchar(k)), collapse = ""), k)
      ggsave(filename = file.path(output_folder_plots_current, paste0(name, ".png")),
             plot = p, #+force_panelsizes(rows = unit(6, "cm"),cols = unit(6, "cm")),
             width = 588, height = 196, units = "px", dpi = 110) # larger dpi = larger font
    }

    expect_true(0 == 0)
  })
}
