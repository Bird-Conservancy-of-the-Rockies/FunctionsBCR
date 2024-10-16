RunNimbleParallel <-
  function(model, inits, data, constants, parameters,
           par.ignore.Rht = c(), par.dontign.Rht = c(),
           par.fuzzy.track.Rht = c(), fuzzy.Rht.threshold = 0.05,
           nc = 2, ni = 2000, nb = 0.5, nt = 10, mod.nam = "mod",
           max.samples.saved = 10000, rtrn.model = F, sav.model = T,
           Rht.required = 1.1, neff.required = 100, check.progress = NULL,
           max.tries = NULL) {
    if(nb < 1 & (ni - (ni * nb)) < 100) stop("Increase iterations (ni) or reduce burn-in. Too few samples for calculating Rhat.")
    if(nb >= 1 & (ni - nb) < 100) stop("Increase iterations (ni) or reduce burn-in. Too few samples for calculating Rhat.")
    
    require(nimble)
    require(parallel)
    require(coda)
    require(mcmcOutput)
    set.seed(1)
    cl<-makeCluster(nc, timeout = 5184000)
    model <<- model
    inits <<- inits
    #data <<- data
    data.inp <<- data
    constants <<- constants
    parameters <<- parameters
    ni <<- ni
    nt <<- nt
    clusterExport(cl, c("model", "inits", "data.inp", "constants", "parameters",
                        "ni", "nt"))
    for (j in seq_along(cl)) {
      set.seed(j)
      if(is.function(inits)) {
        init <<- inits()
      } else {
        init <<- inits
      }
      clusterExport(cl[j], "init")
    }
    out1 <- clusterEvalQ(cl, {
      library(nimble)
      library(coda)
      model <- nimbleModel(code = model, name = "model",
                           constants = constants, data = data.inp,
                           inits = init, calculate = FALSE)
      Cmodel <- compileNimble(model)
      modelConf <- configureMCMC(model, thin = nt)
      # Example code for switching out samplers:
      # modelConf$removeSamplers(c("beta0", "betaVec", "delta0", "dev.delta0", "deltaVec"))
      # modelConf$addSampler(c("beta0", "betaVec"), type = "RW_block", control = list(tries = 10))
      # modelConf$addSampler(c("delta0", "deltaVec"), type = "RW_block", control = list(tries = 10))
      # modelConf$addSampler(c("dev.delta0"), type = "RW_block", control = list(tries = 10))
      modelConf$setMonitors(parameters)
      modelMCMC <- buildMCMC(modelConf)
      CmodelMCMC <- compileNimble(modelMCMC, project = model)
      CmodelMCMC$run(ni, reset = FALSE)
      return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
      gc(verbose = FALSE)
    })
    for(chn in 1:nc) { # nc must be > 1
      ind.keep <- c()
      for(p in 1:length(parameters)) ind.keep <-
          c(ind.keep, which(str_detect(dimnames(out1[[chn]])[[2]], parameters[p]))) %>% unique()
      out1[[chn]] <- out1[[chn]][,ind.keep]
    }
    
    ## Check convergence ##
    out2 <- out1
    ni.saved <- nrow(out2[[1]])
    for(chn in 1:nc) { # nc must be > 1
      
      if(nb < 1) {
        nb.real <- (round(ni.saved * nb)+1)
      } else {
        nb.real <- (round(nb/nt)+1)
      }
      out2[[chn]] <- out2[[chn]][nb.real:ni.saved,]
    }
    out.mcmc <- coda::as.mcmc.list(lapply(out2, coda::as.mcmc))
    
    mod <- mcmcOutput(out.mcmc)
    sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
    sumTab <- sumTab %>%
      as_tibble() %>%
      mutate(Parameter = row.names(sumTab)) %>%
      dplyr::select(Parameter, mean:f)
    
    if(length(par.ignore.Rht) > 0) {
      if(length(par.dontign.Rht) == 0) {
        ind.ignore <- which(str_detect_any(sumTab$Parameter, par.ignore.Rht))
      } else {
        ind.ignore <- which(str_detect_any(sumTab$Parameter, par.ignore.Rht) &
                              !str_detect_any(sumTab$Parameter, par.dontign.Rht))
      }
      if(length(ind.ignore) > 0) sumTab <- sumTab %>% slice(-ind.ignore)
    }
    mxRht <- sumTab %>% pull(Rhat) %>% max() # na.rm = T
    mn.neff <- sumTab %>% pull(n.eff) %>% min() # na.rm = T
    if(is.na(mxRht) | is.na(mn.neff)) {
      write.csv(sumTab, "Model_summary.csv")
      stop("Error: One or more parameters is not being sampled. Check data, initial values, etc., and try again. See 'Model_summary.csv' for parameters missing Rhat or n.eff.")
    }
    Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
    if(length(par.fuzzy.track.Rht) > 0) {
      for(p in 1:length(par.fuzzy.track.Rht)) {
        pfuz <- par.fuzzy.track.Rht[p]
        Rht.fuzzy <- c(Rht.fuzzy,
                       sumTab %>%
                         filter(str_sub(Parameter, 1, nchar(pfuz) + 1) ==
                                  str_c(pfuz, "[")) %>%
                         pull(Rhat))
      }
      rm(p, pfuz)
    }
    if(!is.null(check.progress)) mxRht.1 <- mxRht
    
    mcmc.info <- c(nchains = nc, niterations = ni,
                   burnin = ifelse(nb<1, nb*ni, nb), nthin = nt)
    mod <- list(mcmcOutput = mod, summary = sumTab, mcmc.info = mcmc.info)
    if(sav.model) R.utils::saveObject(mod, mod.nam) # If running all in one.
    
    ## If has not converged, continue sampling
    if(round(mxRht, digits = 1) > Rht.required |
       mn.neff < neff.required |
       (length(Rht.fuzzy) > 0 &
       ((sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = T) - 1) >
       ((length(Rht.fuzzy) - 1) * fuzzy.Rht.threshold)))) { # Disregarding the 1 included earlier to avoid an error in round()
      n.runs <- 1
      R.utils::saveObject(out1, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.
    }
    rm(mod, sumTab, out1, out2, out.mcmc)
    gc(verbose = F)
    while(round(mxRht, digits = 1) > Rht.required |
          mn.neff < neff.required |
          (length(Rht.fuzzy) > 0 &
           ((sum(round(Rht.fuzzy, digits = 1) > Rht.required, na.rm = T) - 1) >
            ((length(Rht.fuzzy) - 1) * fuzzy.Rht.threshold)))) { # Disregarding the 1 included earlier to avoid an error in round()
      n.runs <- n.runs + 1
      print(str_c("Run = ", n.runs, ". Max Rhat = ", mxRht, " and min neff = ", mn.neff,
                  ". Keep on chugging."))
      
      out2 <- clusterEvalQ(cl, {
        CmodelMCMC$run(ni, reset = FALSE, resetMV = TRUE) # Resume sampling.
        return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
        gc(verbose = F)
      })
      for(chn in 1:nc) { # nc must be > 1
        ind.keep <- c()
        for(p in 1:length(parameters)) ind.keep <-
            c(ind.keep, which(str_detect(dimnames(out2[[chn]])[[2]], parameters[p]))) %>% unique()
        out2[[chn]] <- out2[[chn]][,ind.keep]
      }
      R.utils::saveObject(out2, str_c(mod.nam, "_chunk", n.runs)) # Save samples from previous run to drive.
      
      if(nb < 1) {  # Anticipated number of samples to save (assuming half discarded as burn-in).
        ni2 <- round(((ni / nt) * n.runs * nc) * (1 - nb))
      } else {
        ni2 <- round(((ni / nt) * n.runs * nc) - (nb / nt * nc))
      }
      if(ni2 > max.samples.saved) {
        nt2 <- round(1 / (max.samples.saved / ni2)) # Set additional thinning so that saved iterations don't exceed (by too much) max.samples.saved (specified by user).
      } else {
        nt2 <- 1
      }
      
      # Reassemble chain from chunks and apply additional thinning.
      out1 <- R.utils::loadObject(str_c(mod.nam, "_chunk", 1))
      ni.saved <- nrow(out1[[1]])
      for(chn in 1:nc) { # nc must be > 1
        out1[[chn]] <- out1[[chn]][seq(2, ni.saved, by = nt2),] # Starting at 2 because first iteration is NA for some reason.
      }
      for(r in 2:n.runs) {
        out.r <- R.utils::loadObject(str_c(mod.nam, "_chunk", r))
        ni.saved <- nrow(out.r[[1]])
        for(chn in 1:nc) {
          out.r[[chn]] <- out.r[[chn]][seq(1, ni.saved, by = nt2),]
          out1[[chn]] <- rbind(out1[[chn]], out.r[[chn]])
        }
      }
      
      # Discard specified proportion of initial samples as burn-in
      out3 <- out1
      ni.saved <- nrow(out3[[1]])
      for(chn in 1:nc) {
        if(nb < 1) {
          nb.real <- (round(ni.saved * nb)+1)
        } else {
          nb.real <- round((nb / (nt * nt2))+1)
        }
        out3[[chn]] <- out3[[chn]][nb.real:ni.saved,]
      }
      out.mcmc.update <- coda::as.mcmc.list(lapply(out3, coda::as.mcmc))
      
      mod <- mcmcOutput(out.mcmc.update)
      sumTab <- summary(mod, MCEpc = F, Rhat = T, n.eff = T, f = T, overlap0 = T, verbose = F)
      sumTab <- sumTab %>%
        as_tibble() %>%
        mutate(Parameter = row.names(sumTab)) %>%
        dplyr::select(Parameter, mean:f)
      if(length(par.ignore.Rht) > 0) {
        if(length(par.dontign.Rht) == 0) {
          ind.ignore <- which(str_detect_any(sumTab$Parameter, par.ignore.Rht))
        } else {
          ind.ignore <- which(str_detect_any(sumTab$Parameter, par.ignore.Rht) &
                                !str_detect_any(sumTab$Parameter, par.dontign.Rht))
        }
        if(length(ind.ignore) > 0) sumTab <- sumTab %>% slice(-ind.ignore)
      }
      mxRht <- sumTab %>% pull(Rhat) %>% max() # na.rm = T
      mn.neff <- sumTab %>% pull(n.eff) %>% min() # na.rm = T
      Rht.fuzzy <- 1 # Putting in at least one value to avoid error later....
      if(length(par.fuzzy.track.Rht) > 0) {
        for(p in 1:length(par.fuzzy.track.Rht)) {
          pfuz <- par.fuzzy.track.Rht[p]
          Rht.fuzzy <- c(Rht.fuzzy,
                         sumTab %>%
                           filter(str_sub(Parameter, 1, nchar(pfuz) + 1) ==
                                    str_c(pfuz, "[")) %>%
                           pull(Rhat))
        }
        rm(p, pfuz)
      }
      gc(verbose = FALSE)
      
      mcmc.info <- c(nchains = nc, niterations = ni * n.runs,
                     burnin = ifelse(nb<1, round(ni.saved * nb), nb),
                     nthin = nt * nt2)
      if(!is.null(check.progress)) {
        if(n.runs == check.progress & mxRht >= mxRht.1 & mxRht > 2) {
          warn.message <- str_c("Insufficient progress towards convergence - Rhat has not decreased after ", check.progress,
                                " runs. Abandoning model. Recommend checking model output and trying a different model.")
          mod <- list(mcmcOutput = mod, summary = sumTab, mcmc.info = mcmc.info,
                      warning = warn.message)
          if(sav.model) R.utils::saveObject(mod, mod.nam)
          print(warn.message)
          break
        }
      }
      if(!is.null(max.tries)) {
        if(n.runs == max.tries) {
          warn.message <- "Maximum number of runs reached. Abandoning model. Recommend checking model output and trying a different model."
          mod <- list(mcmcOutput = mod, summary = sumTab, mcmc.info = mcmc.info,
                      warning = warn.message)
          if(sav.model) R.utils::saveObject(mod, mod.nam)
          print(warn.message)
          break
        }
      }
      mod <- list(mcmcOutput = mod, summary = sumTab, mcmc.info = mcmc.info)
      if(sav.model) R.utils::saveObject(mod, mod.nam) # If running all in one.
      
      rm(mod, sumTab, out1, out2, out3, out.r, out.mcmc.update)
      gc(verbose = F)
    }
    if(exists("n.runs")) for(r in 1:n.runs) file.remove(str_c(mod.nam, "_chunk", r))
    if(rtrn.model) return(mod)
    stopCluster(cl)
  }
