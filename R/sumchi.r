# Copyright 2014-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of sadists.
#
# sadists is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sadists is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with sadists.  If not, see <http://www.gnu.org/licenses/>.

#
# Created: 2015.02.27
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the cumulants of the sumchi
# distribution. 
sumchi_cumuls <- function(wts,df,ncp=0,order.max=3) {
	subkappa <- mapply(function(w,dd,nn) 
										 { (w ^ (1:order.max)) * chi_cumuls(df=dd,ncp=nn,order.max=order.max) },
										 wts,df,ncp,SIMPLIFY=FALSE)
	kappa <- Reduce('+', subkappa)
	return(kappa)
}
sumchi.support <- function(wts,df,ncp) {
	minv <- ifelse(min(wts) < 0,-Inf,0)
	maxv <- ifelse(max(wts) > 0,Inf,0)
	retval <- c(minv,maxv)
	return(retval)
}

# dsumchi, psumchi, qsumchi, rsumchi#FOLDUP
#' @title The sum of (non-central) chis distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the weighted sum of non-central
#' chis.
#'
#' @details
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares, where \eqn{\nu_i}{v_i}
#' are the degrees of freedom, and \eqn{\delta_i}{delta_i} are the
#' non-centrality parameters.  
#' Let \eqn{w_i} be given constants. Suppose
#' \deqn{Y = \sum_i w_i \sqrt{X_i}.}{Y = sum w_i sqrt(X_i).}
#' Then \eqn{Y}{Y} follows a weighted sum of chis distribution. 
#'
#' @usage
#'
#' dsumchi(x, wts, df, ncp=0, order.max=6, log = FALSE)
#'
#' psumchi(q, wts, df, ncp=0, order.max=6, lower.tail = TRUE, log.p = FALSE)
#'
#' qsumchi(p, wts, df, ncp=0, order.max=6, lower.tail = TRUE, log.p = FALSE)
#'
#' rsumchi(n, wts, df, ncp=0)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param wts the vector of weights. 
#' We do \emph{not} recycle this against \code{x,q,p,n}.
#' @param df the vector of degrees of freedom. 
#' We do \emph{not} recycle this against \code{x,q,p,n}.
#' @param ncp the vector of non-centrality parameters.
#' We do \emph{not} recycle this against \code{x,q,p,n}.
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#'
#' @return \code{dsumchi} gives the density, \code{psumchi} gives the 
#' distribution function, \code{qsumchi} gives the quantile function, 
#' and \code{rsumchi} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dsumchi psumchi qsumchi rsumchi
#' @seealso sum of chi-square distribution functions,
#' \code{\link{dsumchisq}, \link{psumchisq}, \link{qsumchisq}, \link{rsumchisq}},
#' @template etc
#' @examples 
#' wts <- c(1,-3,4)
#' df <- c(100,20,10)
#' ncp <- c(0,2,1)
#' rvs <- rsumchi(128, wts, df, ncp)
#' dvs <- dsumchi(rvs, wts, df, ncp)
#' qvs <- psumchi(rvs, wts, df, ncp)
#' pvs <- qsumchi(ppoints(length(rvs)), wts, df, ncp)
#' @rdname dsumchi
#' @name sumchi
#' @export 
dsumchi <- function(x, wts, df, ncp=0, order.max=6, log = FALSE) {
	kappa <- sumchi_cumuls(wts,df,ncp,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,support=sumchi.support(wts,df,ncp),log=log)
	return(retval)
}
#' @export
psumchi <- function(q, wts, df, ncp=0, order.max=6, lower.tail = TRUE, log.p = FALSE) {
	kappa <- sumchi_cumuls(wts,df,ncp,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,support=sumchi.support(wts,df,ncp),
																		 lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qsumchi <- function(p, wts, df, ncp=0, order.max=6, lower.tail = TRUE, log.p = FALSE) {
	kappa <- sumchi_cumuls(wts,df,ncp,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,support=sumchi.support(wts,df,ncp),
															lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rsumchi <- function(n, wts, df, ncp=0) {
	subX <- mapply(function(w,dd,nn) { w * sqrt(rchisq(n,df=dd,ncp=nn)) },
										 wts,df,ncp,SIMPLIFY=FALSE)
	X <- Reduce('+', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
