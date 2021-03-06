% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.supELtest.R
\name{print.supELtest}
\alias{print.supELtest}
\title{Print a supELtest object}
\usage{
\method{print}{supELtest}(x, digits = max(3L, getOption("digits") - 3L), quiet = FALSE, ...)
}
\arguments{
\item{x}{the result of a call to the \code{supELtest} function}

\item{digits}{significant digits to print, the default value is \code{max(3L, getOption("digits") - 3L)}}

\item{quiet}{a logical indicating whether to reduce the amount of output or not, the default value is \code{FALSE}}

\item{...}{for future method}
}
\description{
Print the maximally selected EL statistics and the p-value of the test.
}
\examples{
library(survELtest)
nocrossings(survival::Surv(hepatitis$time, hepatitis$censor)~
    hepatitis$group, sided = 1)

## OUTPUT:
## Call:
## nocrossings(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~ 
##     hepatitis$group, sided = 1)
## 
## Decision = 1

## A decision value of 1 means the case of crossings or alternative orderings among the 
## survival functions is excluded. Thus, we can proceed to the one-sided test.

result = supELtest(survival::Surv(hepatitis$time, hepatitis$censor)~
             hepatitis$group, sided = 1)
print(result)

## OUTPUT:
## Call:
## supELtest(formula = survival::Surv(hepatitis$time, hepatitis$censor) ~  
##     hepatitis$group, sided = 1)
## 
## One-sided maximally selected EL test statistic = 10.36, p = 0.006
}
\seealso{
\code{\link{hepatitis}}, \code{\link{supELtest}}, \code{\link{summary.supELtest}}
}
