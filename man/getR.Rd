% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variableImportance.R
\name{getR}
\alias{getR}
\title{Coefficent of correlation between observed
and predicted values}
\usage{
getR(object, obs = NULL, newdata = NULL, method = "spearman")
}
\arguments{
\item{object}{An object of \code{\link{Modelling-class}}.}

\item{obs}{an optional vector of observed values}

\item{newdata}{an optional data.frame containing all variables
used in the model. If omitted, the fitted values are used}

\item{method}{the method used for correlation. By default, is used "spearman" method.
For more detail see {\code{\link[stats]{cor}}} function.}
}
\value{
the R value
}
\description{
The correlation coefficent is computed
between observed values (if present) or
against prediction obtained from fitted
values.
}
