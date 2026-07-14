\name{Impute_missing_covs_rf}
\alias{Impute_missing_covs_rf}
\title{Impute missing covariate values.}
\description{Imputes missing covariate values using a Random Forest model informed with specified covariates.}
\usage{
Impute_missing_covs_rf(dat, v.fill, v.inform, v.factor = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame or matrix containing covariate values.}
  \item{v.fill}{String naming the variable in 'dat' with missing values to be imputed.}
  \item{v.inform}{String vector naming the covariates in 'dat' to be used for informing a random forest model to impute missing values.}
  \item{v.factor}{Optional string vector naming columns in 'dat' that should be treated as categorical (run through classification-mode random forest rather than regression) even if not already stored as factors, e.g. character or integer-coded covariates. Include 'v.fill' in this vector if the variable being imputed is itself categorical but not already stored as a factor. Columns already stored as factors do not need to be listed unless they are meant to be treated as categorical predictors alongside non-factor columns; listing them is harmless and their original type and level order are preserved. Defaults to \code{NULL}, in which case only columns already stored as factors are treated as categorical.}
}
\value{Data frame or matrix of same dimensions as 'dat' but with missing values imputed. Columns named in 'v.factor' are restored to their original class (and, for columns already stored as factors, their original level order) after imputation.}
\author{
Quresh S. Latif
}
\note{
R packages required: randomForest, dplyr, stringr

Columns named in 'v.factor' must be simple atomic types (numeric, integer, character, or an unordered factor); ordered factors, dates, and similar multi-class-attribute types are rejected with an error, since their original type cannot be reliably restored after conversion to a plain factor.

Columns named in 'v.factor' are converted to factors and used directly in the random forest model. As with any random forest categorical predictor, a column with more than 53 unique values (e.g. a near-continuous or ID-like field) will cause \code{randomForest} to fail; 'v.factor' is best suited to genuinely categorical variables with a modest number of levels.
}
