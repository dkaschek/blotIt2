#' blotIt out to list converter
#'
#' Takes the output of blotIt2 and converts it to a list by reading out the
#' attributes.
#'
#' @param out 'out' file as produced by \code{alignME()}.
#' @param use_factors Logical, if set to \code{TRUE} the name columns are converted to
#' factors. By this, the plotted order of target will be same order as they are
#' in the source data file fad through \link{read.wide} to \link{alignME}.
#'
#' @export
blotIt_out_to_list <- function(
  out,
  use_factors = T
) {
  out_list <- list(
    aligned = attr(out, "outputs")$aligned,
    scaled = attr(out, "outputs")$scaled,
    prediction = attr(out, "outputs")$prediction,
    original = attr(out, "outputs")$original,
    parameter = attr(out, "outputs")$parameter,
    fixed = attr(out, "fixed"),
    latent = attr(out, "latent")
  )
  if (use_factors)  {
    out_list$aligned$name <- factor(out_list$aligned$name, levels = unique(out_list$aligned$name))
    out_list$scaled$name <- factor(out_list$scaled$name, levels = unique(out_list$scaled$name))
    out_list$prediction$name <- factor(out_list$prediction$name, levels = unique(out_list$prediction$name))
    out_list$original$name <- factor(out_list$original$name, levels = unique(out_list$original$name))
  }

  return(out_list)
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
#' @importFrom lemon facet_rep_wrap
plot_blotIt <- function(
  out,
  plot_points = "aligned",
  plot_line = "aligned",
  spline = F,
  scales = "fixed",
  plot_caption = T,
  ncol = NULL,
  ...
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

  fixed <- out_list$fixed
  latent <- out_list$latent

  plot_list <- out_list

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


  plot_list$prediction <- subset(plot_list$prediction, ...)
  plot_list$original <- subset(plot_list$original, ...)


  plot_list_points <- plot_list[[plot_points]]
  plot_list_line <- plot_list[[plot_line]]

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
    "Errorband: ", used_errors[[plot_line]],"\n",
    "\n",
    "Date: ", Sys.Date())


  ## plot
  if (plot_points == "aligned" & plot_line == "aligned"){
    g <- ggplot(data = plot_list_points, aes(x = time, y = value, group = fixed, color = fixed, fill = fixed))
    g <- g + facet_rep_wrap(~name , scales = scales, repeat.tick.labels = 'left')
  } else {
    g <- ggplot(data = plot_list_points, aes(x = time, y = value, group = latent, color = latent, fill = latent))
    g <- g + facet_rep_wrap(~name * fixed , scales = scales, repeat.tick.labels = 'left', ncol = ncol)
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
    g <- g + geom_smooth(data = plot_list_line, se = FALSE)
  } else {
    g <- g + geom_line(data = plot_list_line, size = 1)
    g <- g + geom_ribbon(data = plot_list_line, aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.3, lty = 0)
  }




  g <- g + scale_color_discrete(name = legend.name) + scale_fill_discrete(name = legend.name)
  g <- g + theme_bw(base_size = 12) + theme(legend.position = "top", legend.key = element_blank(), strip.background = element_rect(color=NA, fill=NA),
                                            axis.line.x = element_line(size = 0.3, colour = "black"),  axis.line.y = element_line(size = 0.3, colour = "black"),
                                            panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(),
                                            panel.grid.minor = element_blank(), panel.border = element_blank(),
                                            panel.background = element_blank(), plot.margin = unit(c(0,0.5,0.5,0.5), "cm"))
  g <- g + xlab("\nTime [min]") + ylab("Signal [a.u.]\n")# + ggtitle("Title Plot2")
  if (plot_caption) {
    g <- g + labs(caption=caption_text)
  }


  return(g)

}
