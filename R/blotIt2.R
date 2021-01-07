globalVariables(c("name", "time", "value", "sigma", "fixed", "latent", "error"))


#' Read data from csv or txt file
#'
#' Reads data stored in a wide format, i.e. description
#' columns with time points, experimental conditions, etc. are
#' followed by each one column per measured target.
#'
#' @param file character, the name of the data file.
#' @param description numeric index vector of the colums containing
#' the description.
#' @param time numeric index of length 1, the time column.
#' @param header a logical value indicating whether the file contains
#' the names of the variables as its first line.
#' @param ... further arguments being passed to \link{read.csv}.
#'
#' @return data frame with columns "name", "time", "value" and other
#' columns describing the measurements.
#'
#' @export
#' @importFrom utils read.csv
read.wide <- function (file, description = 1:3, time = 1, header = TRUE, ...) {

  mydata <- read.csv(file, header = header, ...)
  allnames <- colnames(mydata)

  # Translate characters to numeric if necessary
  if (is.character(description)) ndescription <- which(colnames(mydata) %in% description) else ndescription <- description
  if (is.character(time)) ntime <- which(colnames(mydata) == time) else ntime <- time

  # Check availability of description and time
  if (length(ndescription) < length(description))
    warning("Not all columns proposed by argument 'description' are available in file. \nTaking the available ones.")

  if (length(ntime) == 0)
    stop("File did not contain a time column as proposed by 'time' argument.")

  # Distinguish description data from measurement data
  Description <- mydata[, ndescription]
  restLong <- unlist(mydata[, -ndescription])

  # Create output data frame
  newdata <- data.frame(
    Description,
    name = rep(allnames[-ndescription], each = dim(mydata)[1]),
    value = restLong
  )

  # Remove missing items
  newdata <- newdata[!is.nan(newdata$value), ]
  newdata <- newdata[!is.na(newdata$value), ]


  colnames(newdata)[ntime] <- "time"

  return(newdata)

}


rssModel <- function(parlist, data_fit,
                     model.expr, errmodel.expr, constraint.expr,
                     jac.model.expr, jac.errmodel.expr,
                     deriv = TRUE) {

  with(c(as.list(data_fit), parlist), {

    prediction <- var <- rep(1, length(value))

    prediction[1:length(prediction)] <- eval(model.expr)
    res <- prediction - value
    # value <- prediction
    var[1:length(var)] <- eval(errmodel.expr)^2

    constr <- eval(constraint.expr)

    residuals <- c(res/sqrt(var), var, constr)

    dwres.dp <- dvar.dp <- NULL
    if (deriv) {

      jac <- lapply(1:length(parlist), function(k) {
        v.mod <- v.err <- rep(0, length(value))
        v.mod[1:length(value)] <- eval(jac.model.expr[[k]])
        v.err[1:length(value)] <- eval(jac.errmodel.expr[[k]])
        dwres.dp <- v.mod/sqrt(var) - v.err*res/var
        dvar.dp <- v.err*2*sqrt(var)
        list(dwres.dp, dvar.dp)
      })

      dwres.dp <- lapply(jac, function(j) j[[1]])
      dvar.dp <- lapply(jac, function(j) j[[2]])

    }

    attr(residuals, "dwres.dp") <- dwres.dp
    attr(residuals, "dvar.dp") <- dvar.dp

    return(residuals)

  })

}




#' Align time-course data based on an Mixed-Effects alignment model
#'
#' The function deals primarily with time-course data of different
#' targets which have been measured under different experimental
#' conditions and whose measured values might be on a different
#' scale, e.g. because of different amplification. The algorithm
#' determines the different scaling and estimates the time-course on
#' a common scale.
#'
#'
#' @param data data frame with obligatory columns "name", "time" and
#' "value". Additionally, \code{data} should have further columns,
#' e.g. characterizing experimental conditions (the fixed effects) and
#' sources of variance between data of the same experimental condition
#' (the latent variables). Alternatively, object of class
#' \code{aligned}.
#' @param model character defining the model by which the values in
#' \code{data} can be described, e.g. "ys/sj"
#' @param errmodel character defining a model for the standard
#' deviation of a value, e.g. "sigma0 + value * sigmaR". This model
#' can contain parameters, e.g. "sigma0" and "sigmaR", or numeric
#' variables from \code{data}, e.g. "value" or "time".
#' @param fixed two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameters contained in \code{model}, e.g. "ys", and "name1, ..."
#' refers to variables in \code{data}, e.g. "Condition".
#' The parameters "par1, ..." are determined specific to the levels of
#' "name1, ...".
#' @param latent two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameters contained in \code{model}, e.g. "sj", and "name1, ..."
#' refers to variables in \code{data}, e.g. "Experiment".
#' @param error two-sided formula of the form
#' \code{par1+par2+... ~ name1+name2+...} where "par1, par2, ..." are
#' parameter containd in \code{error}, e.g. "sigma0" and "sigmaR", and
#' "name1, ..." refers to variables in \code{data}. If the same values
#' of "par1, ..." should be assumed for all data, "name1" can be "1".
#' @param log logical indicating whether all parameters are fitted on
#' log-scale.
#' @param normalize logical indicating whether the fixed effect
#' parameter should be normalized to unit mean.
#' @param verbose logical, print out information about each fit
#' @param normalize_input logical, if TRUE the input will be normalized before
#' scaling. see \link{splitData}.
#' @details Alignment of time-course data is achieved by an alignment
#' model which explains the observed data by a function mixing
#' fixed effects, usually parameters reflecting the "underlying"
#' time-course, and latent variables, e.g. scaling parameters taking
#' account for effects like different amplification or loading, etc.
#' Depending on the measurement technique, the data has constant
#' relative error, or constant absolute error or even a combination
#' of those. This error is described by an error function. The error
#' parameters are usually global, i.e. the same parameter values are
#' assumed for all data points.
#'
#' @return Object of class \code{aligned}, i.e. a data frame of the
#' alignment result containing an attribute "outputs":
#'  a list of data frames
#' \describe{
#' \item{prediction}{original data with value and sigma replaced by
#'                   the predicted values and sigmas}
#' \item{scaled}{original data with the values transformed according
#'               to the inverse model, i.e. \code{model} solved for
#'               the first parameter in \code{fixed}, e.g. "ys".
#'               Sigma values are computed by error propagation
#'               from the inverse model equation.}
#' \item{aligned}{the reduced data with the fixed effects and their
#'                uncertainty, only. The result of the alignment
#'                algorithm.}
#' \item{original}{the original data}
#' \item{parameter}{original data augmented by parameter columns.
#'                  Parameters in each row correspond to the levels of
#'                  fixed, latent or error as passed to alignME().
#'                  Used for initialization or parameter values when
#'                  refitting with modified model.}
#' }
#'
#' The estimated parameters are returned by the attribute "parameters".
#' @example inst/examples/example_alignME.R
#' @seealso \link{read.wide} to read data in a wide column format and
#' get it in the right format for \code{alignME()}.
#' \link{plot1}, \link{plot2}, \link{plot3}, \link{plot4} to plot the result of \code{alignME}.
#' @export
#' @import trust
#' @importFrom rootSolve multiroot
#' @importFrom stats D
#' @importFrom utils getParseData
alignME <- function(data, model = "ys/sj", errmodel = "value*sigmaR",
                      fixed = ys~Condition, latent  = sj~Experiment, error = sigmaR~1,
                      log = TRUE, normalize = TRUE, reduce = FALSE, verbose = FALSE) {

  pdata <- NULL
  if (inherits(data, "aligned")) {
    pdata <- attr(data, "outputs")$parameter
    data <- attr(data, "outputs")$original
  }

  # Stop if formulas have the wrong specification
  if (length(as.character(fixed)) < 3) stop("Left and right-hand side of formula 'fixed' is needed")
  if (length(as.character(latent)) < 3) stop("Left and right-hand side of formula 'latent' is needed")
  if (length(as.character(error)) < 3) stop("Left and right-hand side of formula 'err' is needed")

  # Get fixed and latent effects
  fix <- union(c("name", "time"), getSymbols(as.character(fixed)[3]))
  ran <- union("name", getSymbols(as.character(latent)[3]))
  err <- union("name", getSymbols(as.character(error)[3]))


  # Add intercepts
  if (attr(terms(fixed), "intercept") & length(fix) == 2) fix <- c(fix, "1")
  if (attr(terms(latent), "intercept") & length(ran) == 1) ran <- c(ran, "1")
  if (attr(terms(error), "intercept") & length(err) == 1) err <- c(err, "1")


  # Determine to which class parameters belong
  fixedpars <- getSymbols(as.character(fixed)[2])
  latentpars <- getSymbols(as.character(latent)[2])
  errorpars <- getSymbols(as.character(error)[2])

  # Split mydata in independent blocks by levels of name
  data.list <- splitData(data, fixed, latent, normalize_input = normalize_input,
                         log = log)
  ndata <- length(data.list)

  # Targets, parameters and constraints
  targets <- make.unique(
    sapply(data.list, function(d) as.character(d$name)[1]),
    sep = "_"
  )
  for (n in 1:length(targets)) data.list[[n]]$name <- targets[n]
  parameters <- getSymbols(c(model, errmodel), exclude = colnames(data))
  if (normalize) {
    # work with 0.1% uncertainty on normalization constant
    constraint <- paste("1e3*(mean(", fixedpars[1], ") - 1)")
    cstrength <- 1e3
  } else {
    constraint <- "0"
    cstrength <- 0
  }

  # Determine covariates
  covariates <- union(getSymbols(model, exclude = c(fixedpars, latentpars, errorpars)),
                      getSymbols(errmodel, exclude = c(fixedpars, latentpars, errorpars)))
  cat("Covariates:", paste(covariates, sep = ", "), "\n")


  if (length(setdiff(c(fixedpars, latentpars, errorpars), parameters)) > 0)
    stop("Not all paramters are defined in either arguments 'latent', 'fixed' or 'error'")

  names(parameters)[parameters %in% fixedpars] <- "fixed"
  names(parameters)[parameters %in% latentpars] <- "latent"
  names(parameters)[parameters %in% errorpars] <- "error"

  # Replace data values by model prediction in error model
  errmodel <- replaceSymbols("value", paste0("(", model, ")"), errmodel)

  # Go to log-scale with parameters
  if (log) {
    model <- replaceSymbols(parameters, paste0("exp(", parameters, ")"), model)
    errmodel <- replaceSymbols(parameters, paste0("exp(", parameters, ")"), errmodel)
    constraint <- replaceSymbols(parameters, paste0("exp(", parameters, ")"), constraint)
  }

  # Derivs of the models
  dmodel <- deparse(D(parse(text = model), name = fixedpars[1]))
  jac.model <- lapply(parameters, function(p) deparse(D(parse(text = model), name = p)))
  jac.errmodel <- lapply(parameters, function(p) deparse(D(parse(text = errmodel), name = p)))

  # Normalization constraint
  constraint.expr <- parse(text = constraint)

  # Introduce latent, fixed and error into model and error model
  for (n in 1:length(parameters)) {
    model <- replaceSymbols(parameters[n], paste0(parameters[n], "[", names(parameters)[n], "]"), model)
    errmodel <- replaceSymbols(parameters[n], paste0(parameters[n], "[", names(parameters)[n], "]"), errmodel)
    dmodel <- replaceSymbols(parameters[n], paste0(parameters[n], "[", names(parameters)[n], "]"), dmodel)
    jac.model <- lapply(jac.model, function(myjac) replaceSymbols(parameters[n], paste0(parameters[n], "[", names(parameters)[n], "]"), myjac))
    jac.errmodel <- lapply(jac.errmodel, function(myjac) replaceSymbols(parameters[n], paste0(parameters[n], "[", names(parameters)[n], "]"), myjac))
  }

  cat("Model:        ", model, "\n", sep = "")
  cat("Error Model:  ", errmodel, "\n", sep = "")
  model.expr <- parse(text = model)
  dmodel.expr <- parse(text = dmodel)
  errmodel.expr <- parse(text = errmodel)
  jac.model.expr <- lapply(jac.model, function(myjac) parse(text = myjac))
  jac.errmodel.expr <- lapply(jac.errmodel, function(myjac) parse(text = myjac))



  # Loop over all independent blocks
  #out <- parallel::mclapply(1:ndata, function(i) try({
  out <- lapply(1:ndata, function(i) try({
    cat("Target ", i, "/", ndata, ":", targets[i], "\n", sep = "")

    # Data for an independent block
    data_target <- data.list[[i]]
    data_target$sigma <- NaN
    data_target[["1"]] <- "1"

    para_target <- NULL
    if (!is.null(pdata)) para_target <- pdata[pdata$name %in% data_target$name,]

    data_fit.fixed <- do.call(paste_, data_target[ , fix, drop = FALSE])
    data_fit.latent <- do.call(paste_, data_target[, ran, drop = FALSE])


    # Reduce "technical" replicates, i.e multiple values that occur
    # independently of fixed or latent effects
    if (reduce) {
      if (verbose) cat("Analyzing technical replicates ... ")
      groups <- interaction(data_fit.fixed, data_fit.latent)
      if (any(duplicated(groups))) {
        data_target <- do.call(rbind, lapply(unique(groups), function(g) {
          subdata <- data_target[groups == g, ]
          outdata <- subdata[1, ]
          outdata$value <- mean(subdata$value)
          return(outdata)
        }))
        cat("data points that could not be distinguished by either fixed\nor latent variables have been averaged.\n")

      } else {
        if (verbose) cat("none found.\n")
      }
    }


    # Data for the RSS function
    data_fit <- data.frame(
      data_target[, union(c("name", "time", "value", "sigma"), covariates)],
      fixed = do.call(paste_, data_target[ , fix, drop = FALSE]),
      latent = do.call(paste_, data_target[, ran, drop = FALSE]),
      error = do.call(paste_, data_target[, err, drop = FALSE]),
      stringsAsFactors = FALSE
    )




    fixed.levels <- unique(as.character(data_fit$fixed))
    latent.levels <- unique(as.character(data_fit$latent))
    error.levels <- unique(as.character(data_fit$error))
    all.levels <- unlist(lapply(1:length(parameters), function(k) {
      switch(names(parameters)[k],
             fixed = fixed.levels,
             latent = latent.levels,
             error = error.levels)
    }))

    # Initialize parameters
    parsini <- do.call("c", lapply(1:length(parameters), function(n) {

      if (log) v <- 0 else v <- 1

      l <- switch(
        names(parameters)[n],
        fixed = length(fixed.levels),
        latent = length(latent.levels),
        error = length(error.levels)
      )
      p <- as.character(parameters[n])

      parvalues <- structure(rep(v, l), names = rep(p, l))
      if (!is.null(pdata)) {
        parvalues[1:l] <- para_target[[p]][which(!duplicated(data_fit[[names(parameters[n])]]))]
        if (log) parvalues <- log(parvalues)
      }
      return(parvalues)

    }))
    pars.fixed <- NULL

    # Create mask for derivatives d(fixed)/dp, d(latent)/dp, etc.
    mask <- lapply(1:length(parsini), function(k) {
      effect <- names(parameters[match(names(parsini[k]), parameters)])
      mask.vector <- as.numeric(data_fit[[effect]] == all.levels[k])
      return(mask.vector)
    })


    # Function to translate pars in list of parameters to call rssModel()
    resfn <- function(pars, pars.fixed, deriv = TRUE) {

      pars.all <- c(pars, pars.fixed)

      parlist <- lapply(parameters, function(n) {

        subpar <- pars.all[names(pars.all) == n]
        if (n %in% fixedpars) names(subpar) <- fixed.levels
        if (n %in% latentpars) names(subpar) <- latent.levels
        if (n %in% errorpars) names(subpar) <- error.levels

        return(subpar)

      })
      names(parlist) <- parameters

      c(list(res = rssModel(parlist, data_fit,
                            model.expr, errmodel.expr, constraint.expr,
                            jac.model.expr, jac.errmodel.expr,
                            deriv = deriv)),
        parlist)

    }

    # Objective function as being required by trust
    trustfn <- function(pars, pars.fixed, deriv = TRUE) {

      ndata <- nrow(data_fit)
      res.out <- resfn(pars, pars.fixed, deriv = deriv)$res
      wres <- res.out[1:ndata]
      vars <- res.out[(ndata+1):(2*ndata)]
      constraint <- res.out[2*ndata+1]
      dwres.dp <- attr(res.out, "dwres.dp")
      dvar.dp <- attr(res.out, "dvar.dp")


      bessel <- 1 #(nrow(data_fit) - length(pars))/nrow(data_fit)

      value <- sum(wres^2) +  bessel*sum(log(vars)) + constraint^2
      gradient <- NULL
      hessian <- NULL

      if (deriv) {
        #J.res.out.num <- numDeriv::jacobian(function(...) resfn(...)[[1]], pars, pars.fixed = pars.fixed, method = "simple")

        J.res.out <- do.call(cbind, lapply(1:length(pars), function(k) {

          whichpar <- match(names(pars)[k], parameters)
          J.wres <- dwres.dp[[whichpar]]*mask[[k]]
          J.var <- dvar.dp[[whichpar]]*mask[[k]]
          J.constr <- as.numeric(names(pars[k]) == parameters["fixed"])*cstrength/length(fixed.levels)
          if (log) J.constr <- J.constr * exp(pars[k])


          c(J.wres, J.var, J.constr)

        }))

        J.wres <- J.res.out[1:ndata, , drop = FALSE]
        J.vars <- J.res.out[(ndata+1):(2*ndata), , drop = FALSE]
        J.constr <- J.res.out[2*ndata+1, , drop = FALSE]

        gradient <- as.vector(2*wres%*%J.wres + (bessel/vars)%*%J.vars + 2*constraint*J.constr)
        hessian <- 2*t(rbind(J.wres, J.constr))%*%(rbind(J.wres, J.constr))

      }


      prediction <- wres*sqrt(vars) + data_fit$value
      sigma <- sqrt(vars)

      list(value = value, gradient = gradient, hessian = hessian, prediction = prediction, sigma = sigma)


    }

    # The model function, used for inversion by multiroot.
    myfn <- function(pars, parlist) {

      fixed <- 1:length(pars)
      mylist <- c(list(pars, fixed = fixed), as.list(data_fit), parlist)
      names(mylist)[1] <- fixedpars[1]

      values <- with(mylist, eval(model.expr) - data_fit$value)

      return(values)

    }
    jac.myfn <- function(pars, parlist) {

      fixed <- 1:length(pars)
      mylist <- c(list(pars, fixed = fixed), as.list(data_fit), parlist)
      names(mylist)[1] <- fixedpars[1]

      dvalues <- with(mylist, eval(dmodel.expr))

      return(dvalues)

    }



    if (verbose) cat("Starting fit\n")

    # Estimate all parameters
    # time.out <- system.time(myfit <- trust(trustfn, parsini, pars.fixed = pars.fixed, control = list(report.level = -1)))
    # time.out <- system.time(myfit <- trust::trust(trustfn, parsini, pars.fixed = pars.fixed, rinit = 1, rmax = 10))
    # print(time.out)
    myfit <- trust::trust(trustfn, parsini, pars.fixed = pars.fixed, rinit = 1, rmax = 10)

    if (!myfit$converged) warning(paste("Non-converged fit for target", targets[i]))
    res.out <- resfn(myfit$argument, pars.fixed, deriv = FALSE)
    bessel <- sqrt(nrow(data_fit)/(nrow(data_fit) - length(parsini) + normalize))

    # Test identifiability
    sv <- svd(myfit[["hessian"]])[["d"]]
    tol <- sqrt(.Machine$double.eps)
    nonidentifiable <- which(sv < tol*sv[1])
    if (length(nonidentifiable) > 0) {
      warning("Eigenvalue(s) of Hessian below tolerance. Parameter uncertainties might be underestimated.")
    }

    # Collect fitted parameters
    partable <- data.frame(
      name = targets[i],
      level = c(rep(fixed.levels, length(fixedpars)), rep(latent.levels, length(latentpars)), rep(error.levels, length(errorpars))),
      parameter = names(myfit$argument),
      value = myfit$argument,
      sigma = as.numeric(sqrt(diag(2*MASS::ginv(myfit$hessian))))*bessel,
      nll = myfit$value,
      npar = length(myfit$argument) - normalize,
      ndata = nrow(data_target) # + normalize
    )

    if (log) {
      partable$value <- exp(partable$value)
      partable$sigma <- partable$value*partable$sigma
    }


    if (verbose) {
      cat("Estimated parameters on non-log scale:\n")
      print(partable)
      cat("converged:", myfit$converged, ", iterations:", myfit$iterations, "\n")
      cat("-2*LL: ", myfit$value, "on", nrow(data_target) + normalize - length(myfit$argument), "degrees of freedom\n")
    }
    attr(partable, "value") <- myfit$value
    attr(partable, "df") <- nrow(data_target) + normalize - length(myfit$argument)

    # Generate data.frame with model prediction
    data_out_1 <- data_target
    data_out_1$sigma <- myfit$sigma*bessel
    data_out_1$value <- myfit$prediction

    # Generate data.frame with scaled data
    values.scaled <- rep(0, nrow(data_target))
    if (verbose) cat("Inverting model ... ")
    values.scaled <- try(rootSolve::multiroot(myfn, start = values.scaled, jacfunc = jac.myfn, parlist = res.out[-1])$root, silent = TRUE)
    if (verbose) cat("done\n")
    if (inherits(values.scaled, "try-error")) {
      data_out_2 <- NULL
      warning("Rescaling to common scale not possible. Equations not invertible.")
    } else {
      myderiv <- abs(jac.myfn(values.scaled, res.out[-1]))
      sigmas.scaled <- myfit$sigma*bessel/myderiv
      if (log) {
        values.scaled <- exp(values.scaled)
        sigmas.scaled <- values.scaled * sigmas.scaled
      }
      data_out_2 <- data_target
      data_out_2$value <- values.scaled
      data_out_2$sigma <- sigmas.scaled
    }

    # Generate data.frame with alignment result
    nini <- length(fixed.levels)
    data_out_3 <- data_target[!duplicated(data_fit$fixed), intersect(fix, colnames(data_target))]
    data_out_3$value <- res.out[[fixedpars[1]]]
    if (log) data_out_3$value <- exp(data_out_3$value)
    data_out_3$sigma <- as.numeric(sqrt(diag(2*MASS::ginv(myfit$hessian))))[1:nini]*bessel
    if (log) data_out_3$sigma <- data_out_3$value*data_out_3$sigma


    # Augment data_target by fitted parameters
    data_target_aug <- data_target
    for (k in 1:length(parameters)) {
      effect <- names(parameters)[k]
      mylevels <- as.character(data_fit[[effect]])
      index0 <- which(as.character(partable$parameter) == parameters[k])
      index1 <- match(mylevels, as.character(partable$level[index0]))
      index <- index0[index1]
      data_target_aug[[parameters[k]]] <- partable$value[index]
    }

    list(data_out_1, data_out_2, data_out_3, data_target, data_target_aug, partable)


#  }, silent = FALSE), mc.preschedule = FALSE, mc.cores = cores)
  }, silent = FALSE))

  partable <- do.call(rbind, lapply(out, function(o) {
    if (!inherits(o, "try-error")) o[[6]] else NULL
  }))
  attr(partable, "value") <- do.call(sum, lapply(out, function(o) {
    if (!inherits(o, "try-error")) attr(o[[6]], "value") else 0
  }))
  attr(partable, "df") <- do.call(sum, lapply(out, function(o) {
    if (!inherits(o, "try-error")) attr(o[[6]], "df") else 0
  }))



  out <- lapply(1:5, function(i) {
    do.call(rbind, lapply(out, function(o) {
      if (!inherits(o, "try-error")) o[[i]] else NULL
      }))
  })
  names(out) <- c("prediction", "scaled", "aligned", "original", "parameter")

  fixed <- fix[-(1:2)]
  #fixed <- fixed[fixed != "1"]

  latent <- ran[-1]
  #latent <- latent[latent != "1"]
  error <- err[-1]

  myreturn <- out$aligned
  attr(myreturn, "outputs") <- out
  attr(myreturn, "fixed") <- fixed
  attr(myreturn, "latent") <- latent
  attr(myreturn, "error") <- error
  attr(myreturn, "parameters") <- partable
  class(myreturn) <- c("aligned", "data.frame")


  return(myreturn)


}

#' splitData
#'
#' Split data in independent blocks according to fixed and latent
#' variables as being defined for \link{alignME}.
#'
#' @param data data frame with columns "name", "time", "value" and others
#' @param fixed two-sided formula, see \link{alignME}
#' @param latent two-sided formula, see \link{alignME}
#' @param normalize_input logical if set to TRUE, the input data will normalized
#' per latent effect by dividing by the respective mean. Prevents convergence
#' failure on some hardware when the data for different latent effects differ by
#' to many orders of magnitude.
#' @return list of data frames
#' @export
splitData <- function(data, fixed, latent) {

  if (!"1" %in% colnames(data)) {
    data["1"] <- 1
    intercept <- FALSE
  } else {
    intercept <- TRUE
  }
  specific <- c("name", attr(terms(latent), "term.labels"))
  fixed <- c("name", "time", attr(terms(fixed), "term.labels"))

  paste.mod <- function(...) paste(..., sep = "_")

  specific.data <- Reduce(paste.mod, data[specific])
  fixed.data <- Reduce(paste.mod, data[fixed])

  specific.unique <- unique(specific.data)
  fixed.unique <- unique(fixed.data)

  M <- matrix(0, ncol=length(specific.unique), nrow=length(fixed.unique))

  for(i in 1:nrow(data)) {
    myrow <- which(fixed.unique == fixed.data[i])
    mycol <- which(specific.unique == specific.data[i])
    M[myrow, mycol] <- 1
  }

  mylist <- analyzeBlocks(M)


  if (!intercept) data <- data[ , -which(colnames(data) == "1")]
  list_out <- lapply(mylist, function(l) {

    data[fixed.data %in% fixed.unique[l], ]

  })

  # normalize the data
  if (normalize_input) {
    if (log) {
      for (i in seq_len(length(list_out))) {
        list_out[[i]]$value <- list_out[[i]]$value /
          exp(mean(log(list_out[[i]]$value)))
      }
    } else {
      for (i in seq_len(length(list_out))) {
        list_out[[i]]$value <- list_out[[i]]$value / mean(list_out[[i]]$value)
      }
    }
  }

  return(list_out)

}
