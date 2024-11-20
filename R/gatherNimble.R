gatherNimble <- function(read.path, burnin, ni.block, max.samples.saved) {
  require(coda)
  require(mcmcOutput)
  require(stringr)
  
  cNB <- countNimbleBlocks(read.path, burnin, ni.block)
  m <- cNB$m
  nblks <- cNB$nblks
  chns <- unique(m[,"chn"])
  burnin <- cNB$burnin
  burnin.realized <- cNB$burnin.realized
  burnin.needed <- burnin - burnin.realized
  if(burnin.needed < 0) {
    proc$kill_tree()
    stop("Additional burnin needed is negative. countNimbleBlocks burned extra samples and needs to be checked.")
    }
  
  ## Make one matrix for each chain:
  Sys.sleep(5) # To provide time for lingering files to finish writing.
  
  gathr <- lapply(chns, FUN = function (s) {
    blks <- unique(m[m[,"chn"] == s, "blk"])
    lst <- lapply(blks, FUN = function (b) {
      x <- suppressWarnings(system2(command = "lsof", args = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ), stdout = TRUE))
      while(length(x) > 0) {
        Sys.sleep(5)
        x <- suppressWarnings(system2(command = "lsof", args = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ), stdout = TRUE))
        }
      load(file = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ))
      samp
    })
    as.mcmc(do.call(rbind, lst))
  })
  
  ## Make both chains the same length
  nr <- min(sapply(gathr, FUN = nrow))
  gathr.red <- lapply(gathr, FUN = function (mat) {
    as.mcmc(mat[1:nr, ])
  })
  
  ## Apply additional burnin & thinning if needed
  nc <- max(m[,1])
  if(burnin.needed > 0) {
    gathr.red <- lapply(gathr.red, FUN = function (mat) {
      as.mcmc(mat[-c(1:burnin.needed), ])
    })
  }
  chain.length.now <- dim(gathr.red[[1]])[1] / nc
  if(is.null(max.samples.saved)) max.samples.saved <- chain.length.now
  if(max.samples.saved < chain.length.now) {
    ind.sav <- round(seq(1, max.samples.saved, length.out = max.samples.saved))
    additional.thin.rate <- chain.length.now / length(ind.sav)
    gathr.red <- lapply(gathr.red, FUN = function(x) as.mcmc(x[ind.sav,])) %>%
      as.mcmc.list()
  } else {
    gathr.red <- as.mcmc.list(gathr.red)
    additional.thin.rate <- 1
  }
  gc(verbose = FALSE)
  
  out <- mcmcOutput(gathr.red)
  return(mget(c("out", "nblks", "additional.thin.rate")))
}
