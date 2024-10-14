checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c(),
                        spit.summary = FALSE) {
  s <- summary(mcmcOutput, n.eff = TRUE, verbose = FALSE)
  if(length(par.ignore) > 0) {
    if(length(par.dontign) == 0) {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore))
    } else {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore) &
                            !str_detect_any(s$Parameter, par.dontign))
    }
    if(length(ind.ignore) > 0) s <- s %>% slice(-ind.ignore)
  }
  # if(spit.summary) assign(summary.object.name, s, envir = .GlobalEnv)
  result <- max(s$Rhat) <= Rht.required & min(s$n.eff) >= neff.required
  if(spit.summary) {
    return(mget(c("result", "s")))
  } else {
    return(mget(c("result")))
  }
}
