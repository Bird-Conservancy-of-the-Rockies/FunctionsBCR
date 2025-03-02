\name{LoadFilteredWorkspace}
\alias{LoadFilteredWorkspace}
\title{Load R workspace while excluding specified objects.}
\description{Loads R workspace to current environment while specifying objects to exclude.}
\usage{
LoadFilteredWorkspace(workspace, objects.exclude = c(), target.envir)
}
\arguments{
  \item{workspace}{R workspace file address to load.}
  \item{objects.exclude}{Character vector naming objects to exclude from loaded workspace.}
  \item{target.envir}{Envronment to which filtered workspace will be loaded. If current workspace is the target, you will need to assign the current workspace to an object, and then pass that object to 'target.envir'.}
}
\author{
Quresh S. Latif
}
