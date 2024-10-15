RunNimbleParallel <-
  function(model.path, inits, data, constants, parameters,
           par.ignore = c(), par.dontign = c(),
           nc = 2, ni = 2000, nb = 0.5, nt = 10, mod.nam = "mod",
           max.samples.saved = 10000, rtrn.model = F, sav.model = T,
           Rht.required = 1.1, neff.required = 100,
           check.freq = 10, max.tries = NULL, dump.path = "dump") {
    if(!rtrn.model & !sav.model) stop("There is no way for RunNimbleParallel to save output. Set either rtrn.model = TRUE or sav.model = TRUE.")
    if(nb < 1 & (ni - (ni * nb)) < 100) stop("Increase iterations (ni) or reduce burn-in. Too few samples for calculating Rhat.")
    if(nb >= 1 & (ni - nb) < 100) stop("Increase iterations (ni) or reduce burn-in. Too few samples for calculating Rhat.")

    require(nimble)
    require(processx)
    # require(parallel)
    require(coda)
    require(mcmcOutput)
    # Also requires `parallel` package in Linux.
    if(!dir.exists(dump.path)) dir.create(dump.path)
    save(list = c("model.path", "constants", "data", "inits", "parameters", "ni", "nt"), file = str_c(dump.path, "/NimbleObjects.RData"))
    #[Create R script for kicking off nimble run here]. Call it "ModRunScript.R"
    #___________________________________________________________________________#
    writeLines(text = c(
      "require(nimble)",
      "require(stringr)",
      "require(FunctionsBCR)",

      "chn <- commandArgs(trailingOnly = TRUE)[[1]]",
      "dump.path <- commandArgs(trailingOnly = TRUE)[[2]]",
      "path.NimbleWorkspace <- commandArgs(trailingOnly = TRUE)[[3]]",
      
      "load(path.NimbleWorkspace)",
      "source(model.path)",
      "mod <- runNimble(mod.lst = list(model, constants, data, inits, parameters),",
      "n.iter = ni, n.thin = nt, dump.file.path = NULL)",
      "i <- 0",
      "repeat{",
      "i <- i + 1",
      "dump.file.path <- str_c(dump.path, '/mod_chn', chn, '_', i, '.RData')",
      "mod <- runNimble(comp.mcmc = mod$comp.mcmc, n.iter = ni, dump.file.path = dump.file.path)",
      "}"
    ),
    con = str_c(dump.path, "/ModRunScript.R"))
    #___________________________________________________________________________#
    proc <- process$new(command = "parallel",
                        args = c("Rscript", eval(str_c(dump.path, "/ModRunScript.R")),
                                 "{}",
                                 eval(dump.path),
                                 eval(str_c(dump.path, "/NimbleObjects.RData")),
                                 ":::",
                                 1:nc))
    proc
    mod.check.result <- FALSE
    nchecks <- 1
    while(!any(str_detect(list.files(dump.path), "mod_chn"))) {Sys.sleep(10)} # Wait until proc has started writing to file before going on.
    while(ifelse(is.null(max.tries), !mod.check.result, !mod.check.result & nchecks < max.tries)) {
      Sys.sleep(check.freq)
      
      mod.out <- gatherNimble(read.path = dump.path, burnin = nb, ni.block = ni, max.samples.saved = max.samples.saved)
      mod.check <- checkNimble(mod.out$out, Rht.required = Rht.required, neff.required = neff.required,
                               par.ignore = par.ignore, par.dontign = par.dontign,
                               spit.summary = TRUE)
      mod.check.result <- mod.check$result
      nblks <- mod.out$nblks
      thin.additional <- mod.out$additional.thin.rate
      mcmc.info <- c(nchains = nc, niterations = ni*nblks,
                     burnin = ifelse(nb<1, nb*ni*nblks, nb),
                     nthin = nt*thin.additional)
      sumTab.ignore <- mod.check$s
      if(length(par.ignore) > 0) {
        if(length(par.dontign) > 0) {
          sumTab.ignore <- sumTab.ignore %>%
            filter(!str_detect_any(Parameter, par.ignore) | str_detect_any(Parameter, par.dontign))
        } else {
          sumTab.ignore <- sumTab.ignore %>%
            filter(!str_detect_any(Parameter, par.ignore))
        }
      } 
      if(any(is.na(sumTab.ignore$Rhat)) | any(is.na(sumTab.ignore$n.eff))) {
        proc$kill_tree()
        write.csv(sumTab.ignore, str_c("Model_summary_PID",proc$get_pid(),".csv"))
        stop(str_c("Error: One or more parameters is not being sampled.",
                   " Check data, initial values, etc., and try again.",
                   " See 'Model_summary_PID",proc$get_pid(),
                   ".csv' for parameters missing Rhat or n.eff."))
      }
      sumTab <- mod.check$s
      mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info)
      if(sav.model) R.utils::saveObject(mod, mod.nam)
      if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
      if(!mod.check.result) {
        print(str_c("At check = ", nchecks, ". Max Rhat = ", max(sumTab$Rhat),
                    " and min neff = ", min(sumTab$n.eff)))
      } else {
        print(str_c("Model complete at check = ", nchecks, ". Max Rhat = ", max(sumTab$Rhat),
                    " and min neff = ", min(sumTab$n.eff)))
      }
      nchecks <- nchecks + 1
    }
    proc$kill_tree()
    if(!mod.check.result) {
      warn.message <- str_c("Rhat did not decrease after ", nchecks,
                            " checks. Model abandoned before reaching convergence targets.")
      mod <- list(mcmcOutput = mod.out$out, summary = sumTab, mcmc.info = mcmc.info,
                  warning = warn.message)
      if(sav.model) R.utils::saveObject(mod, mod.nam)
      if(rtrn.model) assign("mod", mod.nam, envir = .GlobalEnv)
    }
    unlink(dump.path, recursive = TRUE)
    gc(verbose = FALSE)
  }
