mcmcOutputSubset <- function(mcmcOutput, par.summarize = c(), par.ignore = c()) {
  ## here, "mcmcOutput" is the incoming mcmcOutput object
  ## "mcmcOutput2" is the outgoing (subsetted) mcmcOutput object
  
  library(mcmcOutput)
  library(coda)
  
  nc <- attr(mcmcOutput, "nChains")
  chn.iter <- nrow(mcmcOutput)/nc
  
  m <- chn.iter * matrix(c(0, rep(1:(nc-1), each = 2), nc), nrow = 2) + c(1,0)
  
  ## need to replace actual periods with escaped periods:
  dot.escape <- function (string) {
    gsub(x = string, pattern = "\\.", replace = "\\\\.") 
  }
  
  if(is.null(par.summarize) & is.null(par.ignore))
    stop("One or both of 'par.summarize' and 'par.ignore' need to be specified.")
  
  if(any(par.summarize %in% par.ignore))
    warning("One or more element of 'par.summarize' is in 'par.ignore'. All parameters in 'par.summarize' will be summarized even if they appear in 'par.ignore'.")
  
  if(!is.null(par.summarize)) {
    want.cols <- sapply(par.summarize, FUN = function (s) {
      base <- dot.escape(s)
      mult <- paste0(base, "\\[")  ## maybe a vector or array
      singl <- paste0(base, "$")   ## maybe a scalar
      c(grep(x = colnames(mcmcOutput), pattern = mult, value = TRUE),
        grep(x = colnames(mcmcOutput), pattern = singl, value = TRUE))
    })
    want.cols <- unlist(want.cols)
  }
  
  if(!is.null(par.ignore)) {
    noWant <- sapply(par.ignore, FUN = function (s) {
      base <- dot.escape(s)
      mult <- paste0(base, "\\[")  ## maybe a vector or array
      singl <- paste0(base, "$")   ## maybe a scalar
      c(grep(x = colnames(mcmcOutput), pattern = mult, value = TRUE),
        grep(x = colnames(mcmcOutput), pattern = singl, value = TRUE))
    })
    noWant <- unlist(noWant)
    
    if(!is.null(par.summarize)) {
      want.cols <- c(want.cols, colnames(mcmcOutput)[!colnames(mcmcOutput) %in% noWant])
    } else {
      want.cols <- colnames(mcmcOutput)[!colnames(mcmcOutput) %in% noWant]
    }
  }
  
  if(length(want.cols) == 1) stop("Unfortunately, creation of mcmcOutput object does not work with only one parameter. Either add a parameter to 'par.summarize' or remove a parameter from 'par.ignore'.")
  mcmcOutput.lst <- lapply(as.data.frame(m), FUN = function(v) {
    as.mcmc(mcmcOutput[,][v[1]:v[2], want.cols, drop = FALSE])
  })
  
  names(mcmcOutput.lst) <- NULL
  
  mcmcOutput2 <- mcmcOutput(as.mcmc.list(mcmcOutput.lst))
  stopifnot(all(dimnames(mcmcOutput2)[[2]] == want.cols))
  
  return(mcmcOutput2)
}
