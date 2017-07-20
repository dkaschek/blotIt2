


getSymbols <- function(char, exclude = NULL) {

  char <- char[char!="0"]
  out <- parse(text=char, keep.source = TRUE)
  out <- utils::getParseData(out)
  names <- unique(out$text[out$token == "SYMBOL"])
  if(!is.null(exclude)) names <- names[!names%in%exclude]
  return(names)

}

replaceSymbols <- function(what, by, x) {

  xOrig <- x
  is.not.zero <- which(x!="0")
  x <- x[is.not.zero]

  mynames <- names(x)

  x.parsed <- parse(text = x, keep.source = TRUE)
  data <- utils::getParseData(x.parsed)

  by <- rep(by, length.out=length(what))
  names(by) <- what

  data$text[data$text%in%what] <- by[data$text[data$text%in%what]]
  data <- data[data$token!="expr",]



  breaks <- c(0, which(diff(data$line1) == 1), length(data$line1))

  out <- lapply(1:(length(breaks)-1), function(i) {

    paste(data$text[(breaks[i]+1):(breaks[i+1])], collapse="")

  })

  names(out) <- mynames
  out <- unlist(out)

  xOrig[is.not.zero] <- out

  return(xOrig)


}


## Analyze matrix with zero and one entries for connection components
analyzeBlocks <- function(M) {
  out <- which(apply(M, 1, sum)==0)
  if(length(out) > 0) {
    M <- M[-out,]
    cat("matrix contains zero rows which have been eliminated\n")
  }
  n <- dim(M)[1]
  rcomponents <- list()
  ccomponents <- list()

  counter <- 0
  while(length(unlist(rcomponents))< n) {

    counter <- counter + 1

    if(length(unlist(rcomponents)) == 0) {
      w <- 1
    } else {
      mysample <- (1:n)[-unlist(rcomponents)]
      w <- min(mysample)
    }

    repeat {

      v <- unique(rapply(as.list(w), function(i) which(M[i,] == 1)))
      wnew <- unique(rapply(as.list(v), function(j) which(M[,j] == 1)))
      if(length(wnew) == length(w)) break
      w <- wnew
    }
    rcomponents[[counter]] <- w
    ccomponents[[counter]] <- v


  }

  return(rcomponents)

}


paste_ <- function(...) paste(..., sep = "_")
