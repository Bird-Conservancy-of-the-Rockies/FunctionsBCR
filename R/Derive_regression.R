Derive_regression <- function(y, x, propagate_uncertainty = TRUE) {
  if(length(x) != dim(y)[2]) stop("Length of x must equal ncol(y) for 'Derive_regression'.")
  nsamp <- dim(y)[1]
  nx <- length(x)
  intercept <- slope <- numeric(length = nsamp)
  for(i in 1:nsamp) {
    if(propagate_uncertainty) {
      md <- lm(y[i, ] ~ x)
      mn <- coefficients(md)
      vc <- vcov(md)
      cfs <- MASS::mvrnorm(1, mn, vc)
      m <- cfs["x"]
      b <- cfs["(Intercept)"]
    } else {
      m <- sum((x - mean(x)) * (y[i,] - mean(y[i,]))) / sum((x - mean(x)) ^ 2)
      b <- (sum(y[i,]) - m * sum(x)) / length(y[i, ])
    }
    slope[i] <- m
    intercept[i] <- b
  }
  return(list(intercept = intercept,
              slope = slope))
}
