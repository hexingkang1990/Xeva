% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/access_slot_sensitivity.R
\name{sensitivity}
\alias{sensitivity}
\title{Get sensitivity for an Xeva object}
\usage{
sensitivity(object, type = c("model", "batch"),
  sensitivity.measure = NULL)
}
\arguments{
\item{object}{The \code{Xeva} dataset.}

\item{type}{Sensitivity type (either model or batch).}

\item{sensitivity.measure}{Name of the \code{sensitivity.measure}. Default \code{NULL} will return all sensitivity measures.}
}
\value{
A \code{data.frame} with model or batch ID and sensitivity values.
}
\description{
Given a Xeva object, it will return a \code{data.frame} detailing sensitivity information.
}
\examples{
data(brca)
head(sensitivity(brca, type="batch"))
head(sensitivity(brca, type="model"))
}
