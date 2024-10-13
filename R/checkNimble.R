checkNimble <- function(mcmcOutput, Rht.required = 1.1, neff.required = 100,
                        par.ignore = c(), par.dontign = c()) {
  s <- summary(mcmcOutput, n.eff = TRUE)
  if(length(par.ignore) > 0) {
    if(length(par.dontign) == 0) {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore))
    } else {
      ind.ignore <- which(str_detect_any(s$Parameter, par.ignore) &
                            !str_detect_any(s$Parameter, par.dontign))
    }
    if(length(ind.ignore) > 0) s <- s %>% slice(-ind.ignore)
  }
  if(max(s$Rhat) <= Rht.required & min(s$n.eff) >= neff.required) {
    return(TRUE)
  } else (
    return(FALSE)
  )
}