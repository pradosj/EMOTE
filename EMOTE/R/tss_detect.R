



#' Compute TSS p.value for a given position according to a beta-binomial model
#'
#' @param k.pos number of observed positive event in the tested condition
#' @param k.neg number of observed positive event in the backgroud condition
#' @param n.pos total number of event in the tested condition
#' @param n.neg total number of event in the backgroud condition
#' @param p.lb lower bound on the p-value
#' @return p-value of the enrichment test in pos condition as compared to background condition
#' @import VGAM
#' @export
tss.model <- function(k.pos,k.neg,n.pos,n.neg,p.lb=1e-30) {
  p <- 1 - pbetabinom.ab(k.pos, n.pos, 1+k.neg, 1+n.neg-k.neg)
  pmax(p,p.lb)
}




#' Combine P-values of several replicates with Fisher's method
#'
#' @param P a matrix of P value with as many column as replicate
#' @return the numeric vector of combined p-values
#' @export
#' @importFrom stats pchisq
p.combine <- function(P) {
  x2 <- -2*rowSums(log(P))
  pchisq(x2,2*ncol(P),lower.tail = FALSE)
}







