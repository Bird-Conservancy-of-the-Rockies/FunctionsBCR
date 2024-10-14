checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
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
  # if(spit.summary) assign(summary.object.name, s, envir = .GlobalEnv)
  result <- max(s.ignored$Rhat) <= Rht.required & min(s.ignored$n.eff) >= neff.required
  if(spit.summary) {
    return(mget(c("result", "s")))
  } else {
    return(mget(c("result")))
  }
}
