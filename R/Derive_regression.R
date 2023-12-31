Derive_regression <- function(y, x) {
  if(length(x) != dim(y)[2]) stop("Length of x must equal ncol(y) for 'Derive_regression'.")
  nsamp <- dim(y)[1]
  nx <- length(x)
  intercept <- slope <- numeric(length = nsamp)
  for(i in 1:nsamp) {
    md <- lm(y[i, ] ~ x)
    mn <- coefficients(md)
    vc <- vcov(md)
    cfs <- MASS::mvrnorm(1, mn, vc)
    intercept[i] <- cfs["(Intercept)"]
    slope[i] <- cfs["x"]
  }
  return(list(intercept = intercept,
              slope = slope))
}
