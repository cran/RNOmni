% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BAT.R
\name{WaldStat}
\alias{WaldStat}
\title{Basic Association Score Test}
\usage{
WaldStat(y, g, X)
}
\arguments{
\item{y}{Numeric phenotype vector.}

\item{g}{Genotype vector.}

\item{X}{Model matrix of covariates.}
}
\value{
Numeric matrix, with 1 row per SNP, containing these columns:
\itemize{
  \item "score", the score statistic.
  \item "se", its standard error.
  \item "z", the Z statistic.
  \item "p", the p-value. 
}
}
\description{
Basic Association Score Test
}
