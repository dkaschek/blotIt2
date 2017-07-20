
ggplot <- function(...) ggplot2::ggplot(...)+ theme_blotIt2()


theme_blotIt2 <- function(base_size = 8, base_family = "") {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]

  theme_bw(base_size = base_size, base_family = base_family) +
    theme(line = element_line(colour = gray),
          rect = element_rect(fill = "white", colour = NA),
          text = element_text(colour = black),
          axis.ticks = element_line(colour = black),
          legend.key = element_rect(colour = NA),
          panel.border = element_rect(colour = black),
          panel.grid = element_line(colour = gray, size = 0.2),
          strip.background = element_rect(fill = "white", colour = NA))
}

#' Plot 1
#'
#' Plot original data together with the prediction of the alignment
#' model.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @export
#' @import ggplot2
#' @seealso \link{plot2}, \link{plot3}, \link{plot4}
plot1 <- function(out, ..., residual = FALSE) {

  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  out$prediction$fixed <- do.call(paste_, out$prediction[, fixed, drop = FALSE])
  out$prediction$latent <- do.call(paste_, out$prediction[, latent, drop = FALSE])

  out$original$fixed <- do.call(paste_, out$original[, fixed, drop = FALSE])
  out$original$latent <- do.call(paste_, out$original[, latent, drop = FALSE])

  out$prediction <- subset(out$prediction, ...)
  out$original <- subset(out$original, ...)

  legend.name <- paste(latent, collapse = ", ")

  residuals <- data.frame(
    time = out$prediction$time,
    name = out$prediction$name,
    value = (out$prediction$value - out$original$value)/out$prediction$sigma,
    latent = out$prediction$latent,
    fixed = out$prediction$fixed,
    sigma = 1
  )

  if (!residual) {
    P <- ggplot(out$prediction, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap( ~ name*fixed, scales = "free") +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = .3, lty = 0) +
      geom_line() +
      geom_point(data = out$original) +
      scale_color_discrete(name = legend.name) + scale_fill_discrete(name = legend.name)

  }

  if (residual) {
    P <- ggplot(residuals, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap( ~ name*fixed, scales = "free") +
      geom_step()+
      scale_color_discrete(name = legend.name) + scale_fill_discrete(name = legend.name)

  }

  return(P)

}

plot1mod <- function(out, ...) {

  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  out$prediction$fixed <- do.call(paste_, out$prediction[, fixed, drop = FALSE])
  out$prediction$latent <- do.call(paste_, out$prediction[, latent, drop = FALSE])

  out$original$fixed <- do.call(paste_, out$original[, fixed, drop = FALSE])
  out$original$latent <- do.call(paste_, out$original[, latent, drop = FALSE])

  out$prediction <- subset(out$prediction, ...)
  out$original <- subset(out$original, ...)

  legend.name <- paste(latent, collapse = ", ")

  residuals <- data.frame(
    time = out$prediction$time,
    name = out$prediction$name,
    value = (out$prediction$value - out$original$value)/out$prediction$sigma,
    latent = out$prediction$latent,
    fixed = out$prediction$fixed,
    sigma = 1
  )

 name.fixed <- with(out$prediction, interaction(name, fixed))
  name.fixed.data <- with(out$original, interaction(name, fixed))
  name.fixed.residuals <- with(residuals, interaction(name, fixed))

  P1 <- lapply(1:length(unique(name.fixed)), function(k) {

    subprediction <- out$prediction[name.fixed == unique(name.fixed)[k], ]
    subdata <- out$original[name.fixed.data == unique(name.fixed)[k], ]


    ggplot(subprediction, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap( ~ name*fixed, scales = "free") +
      geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = .3, lty = 0) +
      geom_line() +
      geom_point(data = subdata) +
      scale_color_tableau(name = legend.name) + scale_fill_tableau(name = legend.name) +
      theme(legend.position = "none",
            axis.text.x  = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())

  })

  P2 <- lapply(1:length(unique(name.fixed)), function(k) {

    subresiduals <- residuals[name.fixed.residuals == unique(name.fixed)[k], ]
    ggplot(subresiduals, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
      facet_wrap( ~ name*fixed, scales = "free") +
      geom_hline(yintercept = c(-1, 0, 1), color = "darkgray", lty = c(2, 1, 2)) +
      #geom_point(pch = 45, size = 5)+
      geom_point(size = 1) + geom_line() +
      scale_color_tableau(name = legend.name) + scale_fill_tableau(name = legend.name) +
      ylab("wres") +
      theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_blank())
  })

  morletPlots <- P1
  rawplot <- P2

  g <- ggplotGrob(P1[[1]] + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lwidth <- sum(legend$widths)


  glets <- lapply(morletPlots, ggplotGrob)
  graws <- lapply(rawplot, ggplotGrob)



  rawlet <- function(raw, let, heights=c(3,1)){
    g <- rbind(let, raw)
    panels <- g$layout[grepl("panel", g$layout$name), ]
    g$heights[unique(panels$t)] <- unit(heights, "null")
    g
  }

  combined <- mapply(rawlet, raw = graws, let=glets, SIMPLIFY = FALSE)

  grid.arrange(do.call(function(...) arrangeGrob(..., ncol = NULL), combined), legend, ncol=2, widths = unit.c(unit(1, "npc") - lwidth, lwidth))


}


#' Plot 2
#'
#' Plot the prediction of the alignment model together with the
#' original data transformed to the scale of the predicted
#' time-course.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @seealso \link{plot1}, \link{plot3}, \link{plot4}
#'
#' @export
plot2 <- function(out, ...) {

  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  out$scaled$fixed <- do.call(paste_, out$scaled[, fixed, drop = FALSE])
  out$scaled$latent <- do.call(paste_, out$scaled[, latent, drop = FALSE])

  out$aligned$fixed <- do.call(paste_, out$aligned[, fixed, drop = FALSE])
  out$aligned$latent <- NA

  out$scaled <- subset(out$scaled, ...)
  out$aligned <- subset(out$aligned, ...)

  legend.name <- paste(latent, collapse = ", ")

  P <- ggplot(out$aligned, aes(x = time, y = value, group = latent, color = latent, fill = latent)) +
    facet_wrap( ~ name*fixed, scales = "free") +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = .3, lty = 0) +
    geom_line() +
    geom_point(data = out$scaled) +
    geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), data = out$scaled, width = 0) +
    scale_color_discrete(name = legend.name) + scale_fill_discrete(name = legend.name)


  return(P)

}

#' Plot 3
#'
#'
#' Plot the prediction of the alignment model grouped by fixed
#' effects.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#' @seealso \link{plot1}, \link{plot2}, \link{plot4}
#' @export
plot3 <- function(out, ...) {

  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  out$aligned$fixed <- do.call(paste_, out$aligned[, fixed, drop = FALSE])
  out$aligned$latent <- NA

  out$aligned <- subset(out$aligned, ...)

  legend.name <- paste(fixed, collapse = ", ")

  P <- ggplot(out$aligned, aes(x = time, y = value, group = fixed, color = fixed, fill = fixed)) +
    facet_wrap( ~ name, scales = "free") +
    geom_ribbon(aes(ymin = value - sigma, ymax = value + sigma), alpha = .3, lty = 0) +
    geom_line() +
    geom_point() +
    scale_color_discrete(name = legend.name) + scale_fill_discrete(name = legend.name)


  return(P)


}

#' Plot 4
#'
#'
#' Plot the prediction of the alignment model grouped by fixed
#' effects with a smoothing spline instead of linear interpolation.
#'
#' @param out output of \link{alignME}.
#' @param ... logical expression used for subsetting the data frames,
#' e.g. \code{name == "pERK1" & time < 60}.
#'
#'
#' @seealso \link{plot1}, \link{plot2}, \link{plot3}
#' @export
plot4 <- function(out, ...) {

  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  out$aligned$fixed <- do.call(paste_, out$aligned[, fixed, drop = FALSE])
  out$aligned$latent <- NA

  out$aligned <- subset(out$aligned, ...)

  legend.name <- paste(fixed, collapse = ", ")

  P <- ggplot(out$aligned, aes(x = time, y = value, group = fixed, color = fixed, fill = fixed)) +
    facet_wrap( ~ name, scales = "free") +
    geom_errorbar(aes(ymin = value - sigma, ymax = value + sigma), width = 0) +
    geom_point() +
    geom_smooth(se = FALSE) +
    scale_color_discrete(name = legend.name) + scale_fill_discrete(name = legend.name)


  return(P)


}

