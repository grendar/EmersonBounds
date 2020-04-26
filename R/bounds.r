#'
#' Bounds on true sensitivity, specificity, PPV, NPV of a new test
#'
#' \code{bounds} computes Emerson et al. (2017) bounds on true sensitivity and specificity of a test method wrt an imperfect reference method, as well as (trivially) implied bounds on PPV, NPV.
#' @param TR concordance table of the test T and the reference R; as in Emerson et al. (2017)
#' @param alphaR sensitivity of R (wrt the gold standard)
#' @param betaR specificity of R (wrt the gold standard)
#' @return list of bounds:
#' \describe{
#'   \item{minaT}{minimal value of sensitivity of the tested method}
#'   \item{maxaT}{maximal value of sensitivity}
#'   \item{minbT}{minimal value of specificity}
#'   \item{maxbT}{maximal value of specificity}
#'   \item{minPPV}{minimal value of PPV}
#'   \item{maxPPV}{maximal value of PPV}
#'   \item{minNPV}{minimal value of NPV}
#'   \item{maxNPV}{maximal value of NPV}
#' }
#' @details The TR table of this form is expected:
#' \tabular{ccc}{
#' T vs R \tab 1 \tab 0 \cr
#' 1 \tab n11 \tab n10 \cr
#' 0 \tab n01 \tab n00 \cr
#' }
#'
#' If the disease prevalence theta, as computed by the formula on p. 6 in Emerson et al. (2017) is not in (0,1), list of NAs is returned.
#'
#' @examples
#' ## TR concordance table must have the same form as in Emerson et al. (2017).
#' ##       T\R 1 0
#' ##       1
#' ##       0
#' TR = c(7, 1, 2, 29)
#' dim(TR) = c(2, 2)
#' TR
#' ## Hence, it happened in 2 cases
#' ## that the test method reported '1' (presence of disease)
#' ## and the reference method reported '0' (absence of disease)
#'
#' alphaR = 0.98
#' betaR = 0.96
#'
#' bounds(TR, alphaR, betaR)
#'
#' @references Emerson S.C., Waikar S.S., Fuentes C., Bonventre J.V., Betensky R.A. (2017). Biomarker validation with an imperfect reference: Issues and bounds, Statistical Methods in Medical Research, 27(10):2933-2945.
#'
#' @export
#'
bounds = function(TR, alphaR, betaR)  {
  #
  # sensT
  aT = TR[1,1]/sum(TR[,1])
  # specT
  bT = TR[2,2]/sum(TR[,2])
  #
  piR = sum(TR[,1])/sum(TR)
  theta = (piR + betaR - 1)/(alphaR + betaR - 1)
  #
  psiR = alphaR*theta/piR
  etaR = (betaR*(1-theta))/(1-piR)
  #
  mintaup = max(0, psiR + aT - 1)
  mintaun = max(0, etaR + bT - 1)
  maxtaup = min(psiR, aT)
  maxtaun = min(etaR, bT)
  #
  #
  if ( (theta > 0) & (theta < 1) ) {
    #
    minaT = (mintaup*piR + (1 - etaR - bT + mintaun)*(1-piR))/theta
    maxaT = (maxtaup*piR + (1 - etaR - bT + maxtaun)*(1-piR))/theta
    minbT = ((1 - psiR - aT + mintaup)*piR + mintaun*(1-piR))/(1-theta)
    maxbT = ((1 - psiR - aT + maxtaup)*piR + maxtaun*(1-piR))/(1-theta)
    minPPV = PPV(alpha=minaT, beta=minbT, theta=theta)
    maxPPV = PPV(alpha=maxaT, beta=maxbT, theta=theta)
    minNPV = NPV(alpha=minaT, beta=minbT, theta=theta)
    maxNPV = NPV(alpha=maxaT, beta=maxbT, theta=theta)
    #
    return(c(minaT, maxaT, minbT, maxbT, minPPV, maxPPV, minNPV, maxNPV))
    #
  }
  else {
    #
    cat('Prevalence is not in (0,1)')
    cat('\n')
    return(c(NA,NA,NA,NA,NA,NA,NA,NA))
    #
  }
  #
  #
}


PPV = function(alpha, beta, theta)  {
  #
  # alpha sensitivity
  # beta specificity
  # theta disease prevalence
  #
  #
  if ( (theta > 0) & (theta < 1) & (alpha > 0) & (alpha < 1) & (beta > 0) & (beta < 1)) {
    #
    PPV = (alpha*theta)/(alpha*theta + (1-beta)*(1-theta))
    #
  }
  else {
    return(NA)
  }
  #
}

NPV = function(alpha, beta, theta)  {
  #
  if ( (theta > 0) & (theta < 1) & (alpha > 0) & (alpha < 1) & (beta > 0) & (beta < 1)) {
    #
    NPV = (beta*(1-theta))/((1-alpha)*theta + beta*(1-theta))
    #
  }
  else {
    return(NA)
  }
  #
}


