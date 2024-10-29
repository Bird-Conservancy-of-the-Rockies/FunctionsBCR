runNimble <- function (mod.lst = NULL, comp.mcmc = NULL, n.iter = 1000,
                       n.thin = 1, dump.file.path = NULL) {
  require(nimble)
  stopifnot(sum(c(is.null(mod.lst), is.null(comp.mcmc))) == 1)
  
  if (!is.null(mod.lst)) {
    nm <- nimbleModel(code = mod.lst[[1]], constants = mod.lst[[2]], data = mod.lst[[3]], inits = mod.lst[[4]], calculate = FALSE)
    cat(paste0("Initialization info:\n", nm$initializeInfo(), "\n"))
    cat(paste0("Calculate check: ", nm$calculate(), "\n"))
    nm.conf <- configureMCMC(nm, monitors = mod.lst[[5]], thin = n.thin)
    nm.mcmc <- buildMCMC(nm.conf)
    nm.C <- compileNimble(nm)
    comp.mcmc <- compileNimble(nm.mcmc, project = nm.C)
    samp <- runMCMC(comp.mcmc, niter = n.iter)
  } else {
    if (!is.null(dump.file.path)) {
      samp <- as.matrix(comp.mcmc$mvSamples)
      save(samp, file = dump.file.path)
      rm(samp)
      comp.mcmc$run(niter = n.iter, reset = FALSE, resetMV = TRUE)
    } else {
      comp.mcmc$run(niter = n.iter, reset = FALSE, resetMV = FALSE)
    }
    samp <- as.matrix(comp.mcmc$mvSamples)
  }
  
  return(mget(c("samp", "comp.mcmc")))
}
