colors_SK <- c("#0147A6","#B1170F","#0b8a00","#ffc801","#02a0c5","#ff7802","#65D413","#ffff00",rep("grey",100))

#' PlotNew 1
#'
#' Plot original data together with the prediction of the alignment
#' model in new style.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plotNew2}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
plotNew1 <- function (out, ..., residual = FALSE)
{
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
  residuals <- data.frame(time = out$prediction$time, name = out$prediction$name,
                          value = (out$prediction$value - out$original$value)/out$prediction$sigma,
                          latent = out$prediction$latent, fixed = out$prediction$fixed,
                          sigma = 1)
  if (!residual) {
    P <- ggplot(out$prediction, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap(~name * fixed, scales = "free") +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.3, lty = 0) +
      geom_line() + geom_point(data = out$original) +
      scale_color_manual(values = colors_SK)+
      scale_fill_manual(values = colors_SK)+
      theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                       axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                       panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                       panel.grid.minor = element_blank(), panel.border = element_blank(),
                                       panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
  }
  if (residual) {
    P <- ggplot(residuals, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap(~name * fixed, scales = "free") +
      geom_step() +
      scale_color_discrete(name = legend.name) +
      scale_fill_discrete(name = legend.name) +
      theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                       axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                       panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                       panel.grid.minor = element_blank(), panel.border = element_blank(),
                                       panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) +
      xlab("\nTime") + ylab("Signal [a.u.]\n") + ggtitle("Plot 1")
  }

  print("Plot 1: Data scaled & Prediction scaled")

  return(P)
}

#' PlotNew 2
#'
#' Plot the prediction of the alignment model together with the
#' original data transformed to the scale of the predicted
#' time-course in new style.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
plotNew2 <-function (out, ...)
{
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
    facet_wrap(~name * fixed, scales = "free") +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.2, lty = 0, colour = "grey") +
    geom_line(color = "black") + geom_point(data = out$scaled) +
    geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), data = out$scaled, width = 0) +
    scale_color_manual(values = colors_SK)+
    scale_fill_manual(values = colors_SK)+
    theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                     axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                     panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                     panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) +
    xlab("\nTime") + ylab("Signal [a.u.]\n") + ggtitle("Plot 2")

  print("Plot 2: Data scaled & Prediction aligned")

  return(P)
}

#' PlotNew 3
#'
#' Plot the prediction of the alignment model grouped by fixed
#' effects in new style.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
plotNew3 <-function (out)
{
  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")
  out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  out$aligned$latent <- NA
  legend.name <- paste(fixed, collapse = ", ")
  P <- ggplot(out$aligned, aes(x = time, y = value, group = fixed, color = fixed, fill = fixed)) +
    facet_wrap(~name, scales = "free") +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.3, lty = 0) +
    geom_line() + geom_point() +
    scale_color_manual(values = colors_SK)+
    scale_fill_manual(values = colors_SK)+
    theme_bw(base_size = 18) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                     axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                     panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                     panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) +
    xlab("\nTime") + ylab("Signal [a.u.]\n") + ggtitle("Plot 3")

  print("Plot 3: Data aligned & Prediction aligned")

  return(P)
}

#' PlotDR 1
#'
#' Plot original dose response data together with the prediction of the alignment
#' model in new style.
#'
#' @param out output of \link{alignME} with columns time, dose, ...
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew2}, \link{plotNew3}, \link{plotDR2}, \link{plotDR3}
plotDR1 <- function (out, ...)
{
  fixed <- attr(out, "fixed")[!attr(out, "fixed") == "dose"]
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  if(!rlang::is_empty(fixed)){
  out$prediction$fixed <- do.call(paste0, out$prediction[, fixed, drop = FALSE])
  out$original$fixed <- do.call(paste0, out$original[, fixed, drop = FALSE])
  }
  out$prediction$latent <- do.call(paste0, out$prediction[, latent, drop = FALSE])
  out$original$latent <- do.call(paste0, out$original[, latent, drop = FALSE])

  out$prediction <- subset(out$prediction, ...)
  out$original <- subset(out$original, ...)

  legend.name <- paste(latent, collapse = ", ")

  out$prediction <- out$prediction %>% filter(!dose==0)
  out$original <- out$original %>% filter(!dose==0)

  print("Dose response 1: Data scaled & Prediction scaled")

  P <- ggplot(out$prediction, aes(x = dose, y = value, color = latent, fill = latent)) +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.1, lty = 0) +
    geom_line(size=1) + geom_point(data = out$original) +
    scale_color_manual(values = colors_SK)+
    scale_fill_manual(values = colors_SK)+
    theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                     axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                     panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                     panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) + annotation_logticks(mid = unit(0.1, "cm"),sides="b", alpha = 0.5) +
    scale_y_continuous(limits=c(0,NA))+
    xlab("\nDose [ng/ml]") + ylab("Signal [a.u.]\n") + ggtitle("Dose Response 1")

  if(rlang::is_empty(fixed)){
    P <- P + facet_wrap(~name, scales = "free")
  } else P <- P + facet_wrap(~name * fixed, scales = "free")

  return(P)
}


#' PlotDR 2
#'
#' Plot the dose response prediction of the alignment model together with the
#' original data transformed to the scale of the predicted
#' time-course in new style.
#'
#' @param out output of \link{alignME} with columns time, dose, ...
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew2}, \link{plotNew3}, \link{plotDR1}, \link{plotDR3}
plotDR2 <- function (out, ...)
{
  fixed <- attr(out, "fixed")[!attr(out, "fixed") == "dose"]
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  if(!rlang::is_empty(fixed)){
  out$scaled$fixed <- do.call(paste0, out$scaled[, fixed, drop = FALSE])
  out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  }

  out$scaled$latent <- do.call(paste0, out$scaled[, latent, drop = FALSE])

  out$scaled <- subset(out$scaled, ...)
  out$aligned <- subset(out$aligned, ...)

  legend.name <- paste(latent, collapse = ", ")

  # remove 0 times as log10(0)=-Inf
  out$aligned <- out$aligned %>% filter(!dose==0)
  out$scaled <- out$scaled %>% filter(!dose==0) %>% mutate(ymin = value-sigma) %>% mutate(ymin = replace(ymin, ymin < 0, 0)) # set ymin < 0 to 0 to be able to plot errorbars with limit 0

  print("Dose response 2: Data scaled & Prediction aligned")

  P <- ggplot(out$aligned, aes(x = dose, y = value, color = latent)) +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.2, lty = 0, colour = "grey") +
    geom_line(color = "black") + geom_point(data = out$scaled) +
    geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), data = out$scaled, width = 0) +
    scale_color_manual(values = colors_SK)+
    scale_fill_manual(values = colors_SK)+
    theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                     axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                     panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                     panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) + annotation_logticks(mid = unit(0.1, "cm"),sides="b", alpha = 0.5) +
    scale_y_continuous(limits=c(0,NA))+
    xlab("\nDose [ng/ml]") + ylab("Signal [a.u.]\n") + ggtitle("Dose Response 2")

  if(rlang::is_empty(fixed)){
    P <- P + facet_wrap(~name, scales = "free")
  } else P <- P + facet_wrap(~name * fixed, scales = "free")

  return(P)
}

#' PlotDR 3
#'
#' Plot the dose response prediction of the alignment model grouped by fixed
#' effects in new style.
#'
#' @param out output of \link{alignME} with columns time, dose, ...
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plotNew1}, \link{plotNew2}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}
plotDR3 <- function (out, ...)
{
  fixed <- attr(out, "fixed")[!attr(out, "fixed") == "dose"]
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  if(!rlang::is_empty(fixed)){
  out$aligned$fixed <- do.call(paste0, out$aligned[, fixed, drop = FALSE])
  }
  legend.name <- paste(fixed, collapse = ", ")

  # remove 0 times as log10(0)=-Inf
  out$aligned <- out$aligned %>% filter(!dose==0)

  print("Dose response 3: Data aligned & Prediction aligned")

  P <- ggplot(out$aligned, aes(x = dose, y = value)) +
    scale_color_manual(values = colors_SK)+
    scale_fill_manual(values = colors_SK)+
    theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                     axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                     panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                     panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm")) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100)) + annotation_logticks(mid = unit(0.1, "cm"),sides="b", alpha = 0.5) +
    scale_y_continuous(limits=c(0,NA))+
    xlab("\nDose [ng/ml]") + ylab("Signal [a.u.]\n") + ggtitle("Dose Response 3")

  if(rlang::is_empty(fixed)){
    P <- P + facet_wrap(~name, scales = "free") +
      geom_line(size=1) + geom_point(data = out$aligned, size = 2) +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.2, lty = 0)
  } else {
    P <- P + facet_wrap(~name, scales = "free") +
      geom_line(size=1, aes(color = fixed)) + geom_point(data = out$aligned, size = 2, aes(color = fixed)) +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma, fill = fixed), alpha = 0.2, lty = 0)
  }

  return(P)
}
