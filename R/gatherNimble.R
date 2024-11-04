gatherNimble <- function(read.path, burnin, ni.block, max.samples.saved) {
  require(coda)
  require(mcmcOutput)
  require(stringr)
  
  ##  This could be made more flexible, but for now use a stereotyped approach:
  fl <- list.files(read.path)
  fl <- fl[which(str_detect(fl, "mod_chn"))]
  
  m <- cbind(chn = str_split(fl, "_", simplify = TRUE)[,2] %>%
               str_sub(4, -1) %>% as.integer,
             blk = str_split(str_split(fl, "_", simplify = TRUE)[,3],
                             "\\.", simplify = TRUE)[,1] %>% as.integer)
  nblks <- min(tapply(m[,2], m[,1], max))
  
  rownames(m) <- fl
  m <- m[order(m[,"chn"],m[,"blk"]),,drop = FALSE]
  
  # set a minimum, to implement burn-in
  if(burnin < 1) burnin <- ni.block * nblks * burnin
  burnin.block <- burnin / ni.block
  m <- m[m[,"blk"] > burnin.block,,drop = FALSE]
  
  chns <- unique(m[,"chn"])
  
  ## Make one matrix for each chain:
  Sys.sleep(5) # To provide time for lingering files to finish writing.
  
  gathr <- lapply(chns, FUN = function (s) {
    blks <- unique(m[m[,"chn"] == s, "blk"])
    lst <- lapply(blks, FUN = function (b) {
      x <- system2(command = "lsof", args = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ), stdout = TRUE)
      # x <- try(suppressWarnings(load(file = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ))))
      while(length(x) > 0) {
        Sys.sleep(5)
        x <- system2(command = "lsof", args = paste0(read.path, "/", rownames(m)[m[,"chn"] == s & m[,"blk"] == b] ), stdout = TRUE)
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
  
  ## Apply additional burnin if needed
  nc <- max(m[,1])
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
  
  out <- mcmcOutput(gathr.red)
  return(mget(c("out", "nblks", "additional.thin.rate")))
}
