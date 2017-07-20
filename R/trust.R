#' @import trustOptim
#' @export
trust <- function(objfun, parinit, control = list(), ...) {

  mypars <- parinit
  objfun.out <- NULL


  iterlim <- control[["maxit"]]
  if (is.null(iterlim)) iterlim <- 100

  fn <- function(x, ...) {
    names(x) <- names(parinit)
    mypars <<- x
    objfun.out <<- objfun(x, ...)
    objfun.out[["value"]]
  }

  gr <- function(x, ...) {
    if (any(x != mypars)) value <- fn(x, ...)
    objfun.out[["gradient"]]
  }

  hs <- function(x, ...) {
    if (any(x != mypars)) value <- fn(x, ...)
    as(objfun.out[["hessian"]], "dgCMatrix")
  }

  out <- trustOptim::trust.optim(parinit, fn, gr, NULL, method = "SR1", control = control, ...)

  myfit <- c(
    objfun.out,
    list(
      argument = mypars,
      iterations = out[["iterations"]],
      converged = out[["iterations"]] < iterlim
    )
  )


  return(myfit)


}
