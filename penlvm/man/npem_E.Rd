% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{npem_E}
\alias{npem_E}
\title{E step}
\usage{
npem_E(lambda, psi, Y)
}
\arguments{
\item{lambda}{loading matrix of size p*k we have now}

\item{psi}{the variance matrix of size p*p of common factor we have now}

\item{Y}{the response matrix of size n*p}
}
\value{
non-penalized expectation (ignore the constant)
}
\description{
E step
}
