gatherNimble <- function(read.path, burnin = 1) {
  require(coda)
  require(mcmcOutput)
  
  patt <- ".*_chn([0-9]{1,2})_([0-9]{1,2}).*"
  
  ##  This could be made more flexible, but for now use a stereotyped approach:
  fl <- list.files(read.path)
  fl <- fl[grep(fl, pattern = patt)]
  
  m <- cbind(chn = as.integer(gsub(fl, pattern = patt, replace = "\\1")),
             blk = as.integer(gsub(fl, pattern = patt, replace = "\\2")))
  rownames(m) <- fl
  m <- m[order(m[,"chn"],m[,"blk"]),]
  
  # set a minimum, to implement burn-in
  m <- m[m[,"blk"] > burnin,]
  
  chns <- unique(m[,"chn"])
  
  ## Make one matrix for each chain:
  gathr <- lapply(chns, FUN = function (s) {
    blks <- unique(m[m[,"chn"] == s, "blk"])
    lst <- lapply(blks, FUN = function (b) {
      load(file = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ))
      samp
    })
    as.mcmc(do.call(rbind, lst))
  })
  
  ## to make mcmc.list, all chains must be same length
  ##  if one is out ahead, its latest samples wwill have to be lopped off:
  nr <- sapply(gathr, FUN = nrow)
  # nc <- unique(sapply(gathr, FUN = ncol))
  # stopifnot(length(nc) == 1)
  
  gathr.red <- as.mcmc.list(lapply(gathr, FUN = function (mat) {
    as.mcmc(mat[1:min(nr), ])
  }))
  
  out <- mcmcOutput(gathr.red)
  return(out)
}
