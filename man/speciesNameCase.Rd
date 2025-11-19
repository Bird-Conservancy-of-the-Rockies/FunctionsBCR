\name{speciesNameCase}
\alias{speciesNameCase}
\title{Format species common names to match those in the Bird Conservancy database}
\description{This function formats the capitalization of species common names to match those in the Bird Conservancy database and the standard AOU formatting scheme whereby the first letter of each word in the name is capitalized, and letters after hyphens in the first word are lower case but after hyphens in the second word are capital.}
\usage{
speciesNameCase(rawNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawNames}{Character value or vector of species names with non-standard capitalization.}
}
\value{Character value or vector of species names with standardized capitalization.}
\author{
Quresh S. Latif
}
