% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolsSvenja.R
\name{plotNew3}
\alias{plotNew3}
\title{PlotNew 3}
\usage{
plotNew3(out, ..., ymin = NULL, trans = 0.2, units = c("", ""))
}
\arguments{
\item{out}{output of \link{alignME}}

\item{...}{logical expression used for subsetting the data frames,
e.g. \code{name == "pERK1" & time < 60}.}

\item{ymin}{value defining ymin for all y-axes}

\item{units}{character vector of length 2 defining units for the x- and y-axes,
e.g. \code{c("[min]","[a.u.]")}}
}
\description{
Plot the prediction of the alignment model grouped by fixed
effects in new style.
}
\seealso{
\link{plotNew1}, \link{plotNew3}, \link{plotDR1}, \link{plotDR2}, \link{plotDR3}
}
\author{
Svenja Kemmer
}
