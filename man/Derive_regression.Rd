\name{Derive_regression}
\alias{Derive_regression}
\title{Derive regression for MCMC-sampled parameters}
\description{Applies linear model to series of parameters to derive their relationship with specified covariate. Typical use will be to derive trend for yearly population (e.g., occupancy or abundance) estimates.}
\usage{
Derive_regression(y, x)
}
\arguments{
  \item{y}{Matrix where rows are MCMC samples and columns are parameters to be regressed against x.}
  \item{x}{Vector of covariate values against which y is regressed (should be length = dim(y)[2]).}
}
\value{A list with two elements containing posterior samples of 1) the regression intercept and 2) the regression slope.}
\author{
Quresh Latif
}
\note{R packages required: MASS}

