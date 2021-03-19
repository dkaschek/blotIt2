colors_SK <- c("#0147A6", "#B1170F", "#0b8a00", "#ffc801", "#02a0c5", "#ff7802", "#65D413", "#ffff00", rep("grey", 100))
colors_dMod <- c("#000000", "#C5000B", "#0084D1", "#579D1C", "#FF950E", "#4B1F6F", "#CC79A7","#006400", "#F0E442", "#8B4513", rep("gray", 100))

#' PlotNew 1
#'
#' Plot original data together with the prediction of the alignment
#' model in new style.
#'
#' @param out output of \link{alignME}
#' @param ymin value defining ymin for all y-axes
#' @param units character vector of length 2 defining units for the x- and y-axes,
#' e.g. \code{c("[min]","[a.u.]")}
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}
#'
#' @export
#' @author Svenja Kemmer
#' @import ggplot2
#' @seealso \link{plotNew2}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
plotNew1 <- function(out, ..., ymin = NULL, units = c("", ""), trans = 0.2, residual = FALSE) {
  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")
  out$prediction$fixed <- do.call(paste0, out$prediction[, fixed, drop = FALSE])
  out$prediction$latent <- do.call(paste0, out$prediction[, latent, drop = FALSE])
  out$original$fixed <- do.call(paste0, out$original[, fixed, drop = FALSE])
  out$original$latent <- do.call(paste0, out$original[, latent, drop = FALSE])
  out$prediction <- subset(out$prediction, ...)
  out$original <- subset(out$original, ...)
  legend.name <- paste(latent, collapse = ", ")
  residuals <- data.frame(
    time = out$prediction$time, name = out$prediction$name,
    value = (out$prediction$value - out$original$value) / out$prediction$sigma,
    latent = out$prediction$latent, fixed = out$prediction$fixed,
    sigma = 1
  )
  if (!residual) {
    P <- ggplot(out$prediction, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap(~ name * fixed, scales = "free") +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = trans, lty = 0) +
      geom_line() +
      geom_point(data = out$original) +
      scale_color_manual(values = colors_SK, name = legend.name) +
      scale_fill_manual(values = colors_SK, name = legend.name) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
        axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
      )
  }
  if (residual) {
    P <- ggplot(residuals, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap(~ name * fixed, scales = "free") +
      geom_step() +
      scale_color_discrete(name = legend.name) +
      scale_fill_discrete(name = legend.name) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
        axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
      ) +
      xlab(paste0("Time ", units[1])) +
      ylab(paste0("Signal ", units[2])) +
      ggtitle("Plot 1")
  }

  if (is.numeric(ymin)) P <- P + scale_y_continuous(limits = c(ymin, NA))

  print("Plot 1: Raw data and prediction on original scale")

  return(P)
}

#' PlotNew 2
#'
#' Plot the prediction of the alignment model together with the
#' original data transformed to the scale of the predicted
#' time-course in new style.
#'
#' @param out output of \link{alignME}
#' @param ymin value defining ymin for all y-axes
#' @param units character vector of length 2 defining units for the x- and y-axes,
#' e.g. \code{c("[min]","[a.u.]")}
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @author Svenja Kemmer
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
plotNew2 <- function(out, ..., ymin = NULL, units = c("", "")) {
  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")
  out$scaled$fixed <- do.call(paste0, out$scaled[, fixed, drop = FALSE])
  out$scaled$latent <- do.call(paste0, out$scaled[, latent, drop = FALSE])
  out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  out$aligned$latent <- NA
  out$scaled <- subset(out$scaled, ...)
  out$aligned <- subset(out$aligned, ...)
  legend.name <- paste(latent, collapse = ", ")
  P <- ggplot(out$aligned, aes(x = time, y = value, group = latent, color = latent)) +
    facet_wrap(~ name * fixed, scales = "free") +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.2, lty = 0, colour = "grey") +
    geom_line(color = "black") +
    geom_point(data = out$scaled) +
    geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), data = out$scaled, width = 0) +
    scale_color_manual(values = colors_SK, name = legend.name) +
    scale_fill_manual(values = colors_SK, name = legend.name) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(), panel.border = element_blank(),
      panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    ) +
    xlab(paste0("Time ", units[1])) +
    ylab(paste0("Signal ", units[2])) +
    ggtitle("Plot 2")

  if (is.numeric(ymin)) P <- P + scale_y_continuous(limits = c(ymin, NA))

  print("Plot 2: Data scaled & Prediction aligned")

  return(P)
}

#' PlotNew 3
#'
#' Plot the prediction of the alignment model grouped by fixed
#' effects in new style.
#'
#' @param out output of \link{alignME}
#' @param ymin value defining ymin for all y-axes
#' @param units character vector of length 2 defining units for the x- and y-axes,
#' e.g. \code{c("[min]","[a.u.]")}
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @author Svenja Kemmer
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
plotNew3 <- function(out, ..., ymin = NULL, trans = 0.2, units = c("", "")) {
  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")
  out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  out$aligned$latent <- NA
  out$aligned <- subset(out$aligned, ...)
  legend.name <- paste(fixed, collapse = ", ")
  P <- ggplot(out$aligned, aes(x = time, y = value, group = fixed, color = fixed, fill = fixed)) +
    facet_wrap(~name, scales = "free") +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = trans, lty = 0) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = colors_SK, name = legend.name) +
    scale_fill_manual(values = colors_SK, name = legend.name) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(), panel.border = element_blank(),
      panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    ) +
    xlab(paste0("Time ", units[1])) +
    ylab(paste0("Signal ", units[2])) +
    ggtitle("Plot 3")

  if (is.numeric(ymin)) P <- P + scale_y_continuous(limits = c(ymin, NA))

  print("Plot 3: Data aligned & Prediction aligned")

  return(P)
}

#' PlotDR 1
#'
#' Plot original dose response data together with the prediction of the alignment
#' model in new style.
#'
#' @param out output of \link{alignME} with obligatory columns time and dose
#' @param ymin value defining ymin for all y-axes
#' @param units character vector of length 2 defining units for the x- and y-axes,
#' e.g. \code{c("[ng/ml]","[a.u.]")}
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @author Svenja Kemmer
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew2}, \link{plotNew3}, \link{plotDR2}, \link{plotDR3}
plotDR1 <- function(out, ..., ymin = NULL, units = c("", "")) {
  fixed <- attr(out, "fixed")[!attr(out, "fixed") == "dose"]
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  if (!rlang::is_empty(fixed)) {
    out$prediction$fixed <- do.call(paste0, out$prediction[, fixed, drop = FALSE])
    out$original$fixed <- do.call(paste0, out$original[, fixed, drop = FALSE])
  }
  out$prediction$latent <- do.call(paste0, out$prediction[, latent, drop = FALSE])
  out$original$latent <- do.call(paste0, out$original[, latent, drop = FALSE])

  out$prediction <- subset(out$prediction, ...)
  out$original <- subset(out$original, ...)
  legend.name <- paste(latent, collapse = ", ")

  out$prediction <- out$prediction %>% filter(!dose == 0)
  out$original <- out$original %>% filter(!dose == 0)

  P <- ggplot(out$prediction, aes(x = dose, y = value, color = latent, fill = latent)) +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.1, lty = 0) +
    geom_line(size = 1) +
    geom_point(data = out$original) +
    scale_color_manual(values = colors_SK, name = legend.name) +
    scale_fill_manual(values = colors_SK, name = legend.name) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(), panel.border = element_blank(),
      panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    ) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) +
    annotation_logticks(mid = unit(0.1, "cm"), sides = "b", alpha = 0.5) +
    xlab(paste0("Dose ", units[1])) +
    ylab(paste0("Signal ", units[2])) +
    ggtitle("Dose Response 1")

  if (is.numeric(ymin)) P <- P + scale_y_continuous(limits = c(ymin, NA))

  if (rlang::is_empty(fixed)) {
    P <- P + facet_wrap(~name, scales = "free")
  } else {
    P <- P + facet_wrap(~ name * fixed, scales = "free")
  }

  print("Dose response 1: Data scaled & Prediction scaled")

  return(P)
}


#' PlotDR 2
#'
#' Plot the dose response prediction of the alignment model together with the
#' original data transformed to the scale of the predicted
#' time-course in new style.
#'
#' @param out output of \link{alignME} with obligatory columns time and dose
#' @param ymin value defining ymin for all y-axes
#' @param units character vector of length 2 defining units for the x- and y-axes,
#' e.g. \code{c("[ng/ml]","[a.u.]")}
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @author Svenja Kemmer
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew2}, \link{plotNew3}, \link{plotDR1}, \link{plotDR3}
plotDR2 <- function(out, ..., ymin = NULL, units = c("", "")) {
  fixed <- attr(out, "fixed")[!attr(out, "fixed") == "dose"]
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  if (!rlang::is_empty(fixed)) {
    out$scaled$fixed <- do.call(paste0, out$scaled[, fixed, drop = FALSE])
    out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  }

  out$scaled$latent <- do.call(paste0, out$scaled[, latent, drop = FALSE])

  out$scaled <- subset(out$scaled, ...)
  out$aligned <- subset(out$aligned, ...)
  legend.name <- paste(latent, collapse = ", ")

  # remove 0 times as log10(0)=-Inf
  out$aligned <- out$aligned %>% filter(!dose == 0)
  out$scaled <- out$scaled %>%
    filter(!dose == 0) %>%
    mutate(ymin = value - sigma) %>%
    mutate(ymin = replace(ymin, ymin < 0, 0)) # set ymin < 0 to 0 to be able to plot errorbars with limit 0

  P <- ggplot(out$aligned, aes(x = dose, y = value, color = latent)) +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.2, lty = 0, colour = "grey") +
    geom_line(color = "black") +
    geom_point(data = out$scaled) +
    geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), data = out$scaled, width = 0) +
    scale_color_manual(values = colors_SK, name = legend.name) +
    scale_fill_manual(values = colors_SK, name = legend.name) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(), panel.border = element_blank(),
      panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    ) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) +
    annotation_logticks(mid = unit(0.1, "cm"), sides = "b", alpha = 0.5) +
    xlab(paste0("Dose ", units[1])) +
    ylab(paste0("Signal ", units[2])) +
    ggtitle("Dose Response 2")

  if (is.numeric(ymin)) P <- P + scale_y_continuous(limits = c(ymin, NA))

  if (rlang::is_empty(fixed)) {
    P <- P + facet_wrap(~name, scales = "free")
  } else {
    P <- P + facet_wrap(~ name * fixed, scales = "free")
  }

  print("Dose response 2: Data scaled & Prediction aligned")

  return(P)
}

#' PlotDR 3
#'
#' Plot the dose response prediction of the alignment model grouped by fixed
#' effects in new style.
#'
#' @param out output of \link{alignME} with obligatory columns time and dose
#' @param ymin value defining ymin for all y-axes
#' @param units character vector of length 2 defining units for the x- and y-axes,
#' e.g. \code{c("[ng/ml]","[a.u.]")}
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @author Svenja Kemmer
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew2}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}
plotDR3 <- function(out, ..., ymin = NULL, units = c("", "")) {
  fixed <- attr(out, "fixed")[!attr(out, "fixed") == "dose"]
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  if (!rlang::is_empty(fixed)) {
    out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  }
  out$aligned <- subset(out$aligned, ...)
  legend.name <- paste(fixed, collapse = ", ")

  # remove 0 times as log10(0)=-Inf
  out$aligned <- out$aligned %>% filter(!dose == 0)

  P <- ggplot(out$aligned, aes(x = dose, y = value)) +
    scale_color_manual(values = colors_SK, name = legend.name) +
    scale_fill_manual(values = colors_SK, name = legend.name) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color = NA, fill = NA),
      axis.line.x = element_line(size = 0.3, colour = "black"), axis.line.y = element_line(size = 0.3, colour = "black"),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(), panel.border = element_blank(),
      panel.background = element_blank(), plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    ) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) +
    annotation_logticks(mid = unit(0.1, "cm"), sides = "b", alpha = 0.5) +
    xlab(paste0("Dose ", units[1])) +
    ylab(paste0("Signal ", units[2])) +
    ggtitle("Dose Response 3")

  if (is.numeric(ymin)) P <- P + scale_y_continuous(limits = c(ymin, NA))

  if (rlang::is_empty(fixed)) {
    P <- P + facet_wrap(~name, scales = "free") +
      geom_line(size = 1) + geom_point(data = out$aligned, size = 2) +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.2, lty = 0)
  } else {
    P <- P + facet_wrap(~name, scales = "free") +
      geom_line(size = 1, aes(color = fixed)) + geom_point(data = out$aligned, size = 2, aes(color = fixed)) +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma, fill = fixed), alpha = 0.2, lty = 0)
  }

  print("Dose response 3: Data aligned & Prediction aligned")

  return(P)
}


#' Get scaled replicates
#'
#' Extract scaled replicates from the blotIt output.
#'
#' @param out output of \link{alignME}
#'
#' @export
#' @author Svenja Kemmer
getReplicates <- function(out) {
  attr(out, "outputs")$scaled
}

#' All-in-one plot function for blotIt()
#'
#' Takes the output of 'alignME()' If not already in a list it is internally converted
#' by us of \link{blotIt_out_to_list}. Which data will be plotted can then be
#' specefied seperately.
#'
#' @param out 'out' file as produced by 'alignME()'. Either as raw output or converted
#' to a list, e.g. per \link{blotIt_out_to_list}.
#' @param plot_points String to specify which data set should be plotted in form of points
#' with corresponding error bars. It must be one of \code{c("original", "scaled", "prediction", "aligned")}.
#' @param plot_line Same as above but with a line and errorband.
#' @param spline Logical, if set to \code{TRUE}, what is specified as \code{plot_line} will be plotted
#' as a smooth spline instead of straight lines between points.
#' @param scales String passed as \code{scales} argument to \link{facet_wrap}.
#' @param plot_caption Logical, if \code{TRUE}, a caption describing the plotted data is added to the plot.
#' @param ncol Numerical passed as \code{ncol} argument to \link{facet_wrap}.
#' @param ... Logical expression used for subsetting the data frames, e.g. name == "pERK1" & time < 60.
#'
#'
#' \describe{
#' To reproduce the known function \link{plot1}, \link{plot2} and \link{plot3},
#' use:
#' \item{plot1}{\code{plot_points} = 'original', \code{plot_line} = 'prediction'}
#' \item{plot2}{\code{plot_points} = 'scaled', \code{plot_line} = 'aligned'}
#' \item{plot3}{\code{plot_points} = 'aligned', \code{plot_line} = 'aligned'}
#' }
#'
#'
#' @export
#' @author Severin Bang and Svenja Kemmer
plot_blotItS <- function(
  out,
  ...,
  plot_points = "aligned",
  plot_line = "aligned",
  spline = F,
  scales = "free",
  plot_caption = T,
  ncol = NULL,
  mycolors = NULL,
  duplZeros = FALSE,
  myorder = NULL
) {

  if(!plot_points %in% c("original", "scaled", "prediction", "aligned") |
     !plot_line %in% c("original", "scaled", "prediction", "aligned")) {
    stop("\n\t'plot_points' and 'plot_line' must each be one of c('original', 'scaled', 'prediction', 'aligned')\n")
  }

  if (class(out)[1] == "list") {
    cat("Data is in form of list, continuing\n")
    out_list <- out
  } else {
    cat("Data is not yet in form of a list, passing it thorugh\n
        \t'blotIt_out_to_list(out,use_factors = T)'\n")
    out_list <- blotIt_out_to_list(out,use_factors = T)
  }

  # change plotting order from default
  if(!is.null(myorder)){
    if(length(setdiff(levels(out_list[[1]]$name), myorder)) != 0) {
      stop("myorder doesn't contain all protein names.")
    } else {
      out_list$aligned$name <- factor(out_list$aligned$name, levels = myorder)
      out_list$scaled$name <- factor(out_list$scaled$name, levels = myorder)
      out_list$prediction$name <- factor(out_list$prediction$name, levels = myorder)
      out_list$original$name <- factor(out_list$original$name, levels = myorder)
    }
  }


  fixed <- out_list$fixed
  latent <- out_list$latent

  plot_list <- out_list

  # duplicate 0 values for all doses
  if(duplZeros) {
    for(ndat in 1){
      dat <- plot_list[[ndat]]
      subset_zeros <- copy(subset(dat,time == 0))
      mydoses <- setdiff(unique(dat$dose), 0)
      my_zeros_add <- NULL
      for(d in 1:length(mydoses)){
        subset_zeros_d <- copy(subset_zeros)
        subset_zeros_d$dose <- mydoses[d]
        my_zeros_add <- rbind(my_zeros_add, subset_zeros_d)
      }
      dat <- rbind(dat, my_zeros_add)
      plot_list[[ndat]] <- dat
    }
  }

  # add columns containig the respective latent and fixed effects

  # aligned
  plot_list$aligned$fixed <- do.call(
    paste0,
    plot_list$aligned[, fixed, drop = FALSE]
  )
  plot_list$aligned$latent <- NA

  # scaled
  plot_list$scaled$fixed <- do.call(
    paste0,
    out_list$scaled[, fixed, drop = FALSE]
  )
  plot_list$scaled$latent <- do.call(
    paste0,
    out_list$scaled[, latent, drop = FALSE]
  )

  # prediction
  plot_list$prediction$fixed <- do.call(
    paste0,
    out_list$prediction[, fixed, drop = FALSE]
  )
  plot_list$prediction$latent <- do.call(
    paste0,
    out_list$prediction[, latent, drop = FALSE]
  )

  # original
  plot_list$original$fixed <- do.call(
    paste0,
    out_list$original[, fixed, drop = FALSE]
  )
  plot_list$original$latent <- do.call(
    paste0,
    out_list$original[, latent, drop = FALSE]
  )

  plot_list_points <- plot_list[[plot_points]]
  plot_list_line <- plot_list[[plot_line]]

  plot_list_points <- subset(plot_list_points, ...)
  plot_list_line <- subset(plot_list_line, ...)

  legend.name <- paste(latent, collapse = ", ")

  # build Caption
  used_errors <- list(
    aligned = "Fisher Information",
    scaled = "Propergated error model to common scale",
    prediction = "Error model",
    original = "None"
  )

  used_data <- list(
    aligned = "Estimated true values",
    scaled = "Original data scaled to common scale",
    prediction = "Predictions from model evaluation on original scale",
    original = "Original data"
  )

  caption_text <- paste0(
    "Datapoints: ", used_data[[plot_points]], "\n",
    "Errorbars: ", used_errors[[plot_points]], "\n",
    "Line: ", used_data[[plot_line]],"\n",
    if(plot_points != plot_line) paste0("Errorband: ", used_errors[[plot_line]],"\n"),
    "\n",
    "Date: ", Sys.Date())

  # we want to keep the x ticks!
  if(scales == "fixed") scales <- "free_x"

  if(is.null(mycolors)) mycolors <- colors_dMod
  else mycolors <- c(mycolors, rep("gray", 100))

  ## plot
  if (plot_points == "aligned" & plot_line == "aligned"){
    g <- ggplot(data = plot_list_points, aes(x = time, y = value, group = fixed, color = fixed, fill = fixed))
    g <- g + facet_wrap(~name , scales = scales, ncol = ncol)
    g <- g + scale_color_manual("Condition", values = mycolors) + scale_fill_manual("Condition", values = mycolors)
  } else {
    g <- ggplot(data = plot_list_points, aes(x = time, y = value, group = latent, color = latent, fill = latent))
    g <- g + facet_wrap(~name * fixed , scales = scales, ncol = ncol)
    g <- g + scale_color_manual("Scaling", values = mycolors) + scale_fill_manual("Scaling", values = mycolors)
  }
  g <- g + geom_point(data = plot_list_points, size =2)
  g <- g + geom_errorbar(data = plot_list_points, aes(ymin = value - sigma, ymax = value + sigma), size=1, width = 4)


  if (spline) {
    g <- g + geom_errorbar(
      data = plot_list_line,
      aes(
        ymin = value - sigma,
        ymax = value + sigma
      ),
      width = 0)
    g <- g + geom_smooth(data = plot_list_line, se = FALSE, method = 'lm', formula = y ~ poly(x, 3))
  } else {
    if (plot_points == plot_line | plot_line == "prediction"){
      g <- g + geom_line(data = plot_list_line, size = 1)
      if(plot_line == "prediction") g <- g + geom_ribbon(data = plot_list_line, aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.1, lty = 0)
    } else {
      g <- g + geom_line(data = plot_list_line, size = 1, color = "grey")
      g <- g + geom_ribbon(data = plot_list_line, aes(ymin = value - sigma, ymax = value + sigma, fill = "grey", color = "grey"), alpha = 0.3, lty = 0)
    }
  }

  g <- g + theme_bw(base_size = 20) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                            axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                            panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                            panel.grid.minor = element_blank(), panel.border = element_blank(),
                                            panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
  g <- g + xlab("\nTime") + ylab("Signal\n")# + ggtitle("Title Plot2")

  if(plot_points != "original"){
    # scale y-axes (let them start at same minimum determined by smallest value-sigma and end at individual ymax)
    plot_list_points <- as.data.table(plot_list_points)
    blank_data <- plot_list_points[, list(ymax = max(value + sigma), ymin = min(value - sigma)), by = c("name", "fixed", "latent")]
    blank_data[, ":=" (ymin = min(ymin), ymax = ymaximal(ymax)), by = c("name", "fixed", "latent")]
    blank_data <- melt(blank_data, id.vars = c("name", "fixed", "latent"), measure.vars = c("ymax", "ymin"), value.name = "value") %>% .[, ":="(time = 0, variable = NULL)]
    g <- g + geom_blank(data = as.data.frame(blank_data), aes(x = time, y = value))
  }

  if (plot_caption) {
    g <- g + labs(caption=caption_text)
  }


  return(g)

}

ymaximal <- function(x){
  rx <- round(x, 1)
  lower <- floor(x)
  upper <- ceiling(x)
  lowdiff <- rx-lower
  uppdiff <- upper-rx
  if(x > 5){
    if(upper %% 2 == 0){out <- upper} else
      out <- lower + 1
  } else
    if(x < 2){
      if(rx > x){
        if(rx %% 0.2 == 0){out <- rx} else
          out <- rx + 0.1} else
            if(rx %% 0.2 != 0){out <- rx + 0.1} else
              out <- rx + 0.2
    } else
      if(lowdiff < uppdiff){out <- lower + 0.5} else
        out <- upper
  return(out)
}


