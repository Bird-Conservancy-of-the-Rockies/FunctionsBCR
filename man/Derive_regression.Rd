\name{Derive_regression}
\alias{Derive_regression}
\title{Derive regression for MCMC-sampled parameters}
\description{Applies linear model to derive the relationship of a parameter (y) estimated with Bayesian methods for a set of observations (n) with a specified covariate (x). Typical use will be to derive trend for yearly population (e.g., occupancy or abundance) estimates.}
\usage{
Derive_regression(y, x, propagate_uncertainty = TRUE)
}
\arguments{
  \item{y}{Matrix where rows are MCMC samples and columns are observations. Number of columns must be equal in length to x.}
  \item{x}{Vector of covariate values against which y is regressed (should be length = dim(y)[2]).}
  \item{propagate_uncertainty}{Logical indicating whether to draw coefficient values for each sample from a multivariate normal distribution representing uncertainty in the intercept and slope estimates. If TRUE (default), the stats::lm function is applied to fit a general linear model, and intercept and slope estimates are drawn from a multivariate normal defined by the variance-covariate matrix representing model-estimated uncertainty. If FALSE, intercept and slope parameters are calculated using formulae for least squares regression (Neter et al. 1996:19).}
}
\value{A list with two elements containing posterior samples of 1) the regression intercept and 2) the regression slope.}
\author{
Quresh Latif
}
\note{R package MASS required if propagate_uncertainty = TRUE}
\references{Neter, J., M. H. Kutner, C. J. Nachtsheim, and W. Wasserman. 1996. Applied linear statistical models. Times Mirror Higher Education Group, Inc.}

