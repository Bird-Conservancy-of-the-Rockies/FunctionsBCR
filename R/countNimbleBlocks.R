countNimbleBlocks <- function(read.path, burnin, ni.block) {
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
  
  return(mget(c("m", "nblks")))
}
