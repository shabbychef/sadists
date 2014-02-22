#
# Created: 2014.02.21
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' @title A function.
#'
#' @description 
#'
#' One sentence or so that tells you some more.
#'
#' @details
#'
#' Really detailed. \eqn{\zeta}{zeta}.
#'
#' A list:
#' \itemize{
#' \item I use \eqn{n}{n} to stand for blah.
#' \item and so forth....
#' }
#'
#' @usage
#'
#' funcname(x, n, zeta, ...)
#'
#' @param x vector of blah
#' @param n number of blah
#' @param ... arguments passed on to ...
#' @template param-ope
#' @return \code{dsr} gives the density, \code{psr} gives the distribution function,
#' \code{qsr} gives the quantile function, and \code{rsr} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @keywords distribution 
#' @keywords io
#' @keywords plotting
#' @aliases psr qsr rsr
#' @seealso t-distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
#' @note
#' This is a thin wrapper on the t distribution. 
#' @export 
#' @template etc
#' @template sr
#' @examples 
#' rvs <- rsr(128, 253*6, 0, 253)
#' dvs <- dsr(rvs, 253*6, 0, 253)
#' pvs.H0 <- psr(rvs, 253*6, 0, 253)
#' pvs.HA <- psr(rvs, 253*6, 1, 253)
#' \dontrun{
#' plot(ecdf(pvs.H0))
#' plot(ecdf(pvs.HA))
#' }
#' @author Steven E. Pav \email{steven@@cerebellumcapital.com}

#' @references 
#'
#' Lecoutre, Bruno. "Another look at confidence intervals for the noncentral t distribution." 
#' Journal of Modern Applied Statistical Methods 6, no. 1 (2007): 107--116.
#' \url{http://www.univ-rouen.fr/LMRS/Persopage/Lecoutre/telechargements/Lecoutre_Another_look-JMSAM2007_6(1).pdf}
#'
#' Lecoutre, Bruno. "Two useful distributions for Bayesian predictive procedures under normal models."
#' Journal of Statistical Planning and Inference 79  (1999): 93--105. 
