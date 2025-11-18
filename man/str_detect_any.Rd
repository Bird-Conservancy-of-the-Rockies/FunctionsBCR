\name{str_detect_any}
\alias{str_detect_any}
\title{Vector to vector string matching}
\description{Wrapper function of `stringr::str_detect` that detects which elements of a given vector contain any element of a second string vector.}
\usage{
str_detect_any(string, pattern.vec, ignore_case = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{string}{Character vector in which to search for a set of patterns.}
  \item{pattern.vec}{Character vector listing whose elements are patterns to search for in 'string'.}
  \item{ignore_case}{If TRUE, ignore case. Default is FALSE.}
}
\value{Logical vector indicating which elements of 'string' contain one or more of the elements in 'pattern.vec.}
\author{
Quresh S. Latif
}
\note{R packages required: stringr, dplyr}