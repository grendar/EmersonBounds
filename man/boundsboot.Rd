% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boundsboot.r
\name{boundsboot}
\alias{boundsboot}
\title{Statistic, to be used in boot() function of the library boot}
\usage{
boundsboot(data, indices, abR)
}
\arguments{
\item{data}{data frame of results of the T test, R test; in this order}

\item{indices}{needed by the boot() function}

\item{abR}{vector of sensitivity alphaR, specificity betaR of the reference method (wrt the gold standard)}
}
\value{
list of bounds:
\describe{
  \item{minaT}{minimal value of sensitivity of the tested method}
  \item{maxaT}{maximal value of sensitivity}
  \item{minbT}{minimal value of specificity}
  \item{maxbT}{maximal value of specificity}
  \item{minPPV}{minimal value of PPV}
  \item{maxPPV}{maximal value of PPV}
  \item{minNPV}{minimal value of NPV}
  \item{maxNPV}{maximal value of NPV}
}
to be used by the boot() function.
}
\description{
\code{boundsboot} is the statistic, to be used in boot() function of the library boot, for obtaining the bootstrap confidence interval for the lower and upper bound on the true sensitivity, specificity, PPV and NPV of the test T
}
\details{
In some settings, some bootstrap samples may lead disease prevalence theta outside (0,1).
}
\examples{

library(boot)

data('bbb', package = 'EmersonBounds')

alphaR = 0.98
betaR = 0.96

set.seed(321)
boo = boot(data=bbb, statistic=boundsboot, R=1000, abR=c(alphaR, betaR))

plot(boo, index=4)      # bootstrap diagnostic plot for the upper bound on specificity

boot.ci(boo, index=4)   # 95\% boot ci for the upper bound on specificity

## Try also, e.g.,
## alphaR = 0.91
## betaR = 0.88
## and other settings to see that bootstrap samples
## may lead TR that gives prevalence outside (0,1)
##




}
\references{
Emerson S.C., Waikar S.S., Fuentes C., Bonventre J.V., Betensky R.A. (2017). Biomarker validation with an imperfect reference: Issues and bounds, Statistical Methods in Medical Research, 27(10):2933-2945.
}
