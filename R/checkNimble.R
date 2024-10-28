checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
                        par.fuzzy.track = c(), fuzzy.threshold = 0.05,
                        spit.summary = FALSE) {
  s <- summary(mcmcOutput, MCEpc = T, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = FALSE)
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
    if(length(ind.ignore) > 0) s.ignored <- s %>% slice(-ind.ignore)
  }
  result <- max(s.ignored$Rhat) <= Rht.required & min(s.ignored$n.eff) >= neff.required
  if(length(par.fuzzy.track) > 0) {
    Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
    for(p in 1:length(par.fuzzy.track)) {
      pfuz <- par.fuzzy.track[p]
      Rht.fuzzy <- c(Rht.fuzzy,
                     s %>% filter(str_sub(Parameter, 1, nchar(pfuz) + 1) == str_c(pfuz, "[")) %>%
                       pull(Rhat))
    }
    if(((sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = T) - 1) >
        ((length(Rht.fuzzy) - 1) * fuzzy.threshold))) result <- FALSE
  }
  if(spit.summary) {
    return(mget(c("result", "s")))
  } else {
    return(mget(c("result")))
  }
}
