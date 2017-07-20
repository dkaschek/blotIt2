#' @export
parameters <- function(out) {
  attr(out, "parameters")
}

#' @export
"parameters<-" <- function(out, value) {
  attr(out, "parameters") <- value
  return(out)

}

#' @export
outputs <- function(out) {
  attr(out, "outputs")
}

#' @export
"outputs<-" <-  function(out, value) {
  attr(out, "outputs") <- value
  return(out)
}

#' @export
subset.aligned <- function(x, subset, ...) {

  e <- substitute(subset)

  outputs <- lapply(outputs(x), function(o) o[eval(e, o, parent.frame()),])
  parameters <- parameters(x)[eval(e, parameters(x), parent.frame()), ]
  x <- x[eval(e, x, parent.frame()), ]
  outputs(x) <- outputs
  parameters(x) <- parameters
  class(x) <- c("aligned", "data.frame")

  return(x)

}

#' @export
llrtest <- function(H0, H1, check = TRUE) {

  fixed0 <- union(attr(H0, "fixed"), "1")
  fixed1 <- union(attr(H1, "fixed"), "1")
  latent0 <- union(attr(H0, "latent"), "1")
  latent1 <- union(attr(H1, "latent"), "1")
  error0 <- union(attr(H0, "error"), "1")
  error1 <- union(attr(H1, "error"), "1")

  if (check) {
    if (!all(fixed0 %in% fixed1) | !all(latent0 %in% latent1) | !all(error0 %in% error1))
      stop("H0 is not a special case of H1.")

  }


  value0 <- attr(parameters(H0), "value")
  value1 <- attr(parameters(H1), "value")
  df0 <- attr(parameters(H0), "df")
  df1 <- attr(parameters(H1), "df")

  list(
    llr = value0 - value1,
    statistic = paste0("chisquare with ", df0 - df1, " degrees of freedom."),
    p.value = pchisq(value0 - value1, df = df0 - df1, lower.tail = FALSE)
  )


}

#' @export
residuals <- function(out, ...) {

  fixed <- attr(out, "fixed")
  latent <- attr(out, "latent")
  out <- attr(out, "outputs")

  out$prediction$fixed <- do.call(paste_, out$prediction[, fixed, drop = FALSE])
  out$prediction$latent <- do.call(paste_, out$prediction[, latent, drop = FALSE])

  out$original$fixed <- do.call(paste_, out$original[, fixed, drop = FALSE])
  out$original$latent <- do.call(paste_, out$original[, latent, drop = FALSE])

  out$prediction <- subset(out$prediction, ...)
  out$original <- subset(out$original, ...)


  data.frame(
    time = out$prediction$time,
    name = out$prediction$name,
    value = (out$prediction$value - out$original$value)/out$prediction$sigma,
    latent = out$prediction$latent,
    fixed = out$prediction$fixed,
    sigma = 1
  )

}
