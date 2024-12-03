checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
                        par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                        spit.summary = FALSE, mod.nam = "mod") {
  require(mcmcOutput)
  
  if(!is.null(par.ignore)) {
    ind.cols.check <- which(!str_detect_any(colnames(mcmcOutput), par.ignore))
    if(!is.null(par.dontign)) {
      ind.cols.check <- c(ind.cols.check,
                          which(str_detect_any(colnames(mcmcOutput), par.dontign))) %>%
        unique()
      if(!is.null(par.fuzzy.track)) ind.cols.check <- c(ind.cols.check,
                                                        which(str_detect_any(colnames(mcmcOutput), par.fuzzy.track))) %>%
          unique()
    }
    nc <- dim(mcmcOutput[,,])[2]
    mcmcOutput.reduce <- list()
    for(i in 1:nc) mcmcOutput.reduce[[i]] <- as.mcmc(mcmcOutput[,i,ind.cols.check])
    mcmcOutput <- mcmcOutput.reduce %>% as.mcmc.list() %>% mcmcOutput()
    rm(i, mcmcOutput.reduce)
  }
  s <- summary(mcmcOutput, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
  s <- s %>%
    as_tibble() %>%
    mutate(Parameter = row.names(s)) %>%
    dplyr::select(Parameter, mean:f)
  if(length(par.ignore) > 0) {
    if(length(par.dontign) == 0) {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore))
    } else {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore) &
                            !str_detect_any(s$Parameter, par.dontign))
    }
    if(length(ind.ignore) > 0) {
      s.ignored <- s %>% slice(-ind.ignore)
    } else {
      s.ignored <- s
    }
    if(any(is.na(s.ignored$Rhat))) {
      write.csv(s.ignored %>% filter(is.na(Rhat) | Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
      stop(paste0("Parameters missing Rhat. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
    }
    if(any(s.ignored$Rhat %in% c(Inf, -Inf))) {
      write.csv(s.ignored %>% filter(Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_", mod.nam, ".csv"))
      stop(paste0("Parameters with Inf or -Inf Rhats. Check Bad_pars_", mod.nam, ".csv and possibly try alternative initial values or check data."))
    }
    result <- max(s.ignored$Rhat) <= Rht.required & min(s.ignored$n.eff) >= neff.required
  } else {
    result <- max(s$Rhat) <= Rht.required & min(s$n.eff) >= neff.required
  }
  if(length(par.fuzzy.track) > 0) {
    Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
    if(!any(names(s) == "Rhat")) {
      proc$kill_tree()
      write.csv(s, str_c(species, "_sum_at_fail.csv"), row.names = FALSE)
      stop("Stopped model run because Rhat not calculated.")
    }
    for(p in 1:length(par.fuzzy.track)) {
      pfuz <- par.fuzzy.track[p]
      Rht.fuzzy <- c(Rht.fuzzy,
                     s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                       pull(Rhat))
    }
    Rht.fuzzy <- Rht.fuzzy[-1]
    if(any(is.na(Rht.fuzzy))) {
      write.csv(s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                  filter(is.na(Rhat) | Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy parameters missing Rhat. Check Bad_pars_fuzzy_", mod.nam,
                  ".csv and possibly try alternative initial values or check data."))
    }
    if(any(s.ignored$Rhat %in% c(Inf, -Inf))) {
      write.csv(s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                  filter(Rhat %in% c(Inf, -Inf)), paste0("Bad_pars_fuzzy_", mod.nam, ".csv"))
      stop(paste0("Fuzzy parameters with Inf or -Inf Rhats. Check Bad_pars_fuzzy_", mod.nam,
                  ".csv and possibly try alternative initial values or check data."))
    }
    if(((sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = TRUE) + sum(is.na(Rht.fuzzy))) / length(Rht.fuzzy)) >
        (length(Rht.fuzzy) * fuzzy.threshold)) result <- FALSE
  }
  if(spit.summary) {
    return(mget(c("result", "s")))
  } else {
    return(mget(c("result")))
  }
}
