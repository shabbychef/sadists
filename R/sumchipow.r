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

# compute the cumulants of the sumchipow
# distribution. 
sumchipow_cumuls <- function(wts,df,ncp,pow,order.max=3) {
	#nterms <- max(vapply(list(wts,df,ncp),length,0))
	subkappa <- mapply(function(w,dd,nn,pp) 
										 { (w ^ (1:order.max)) * chipow_cumuls(df=dd,ncp=nn,pow=pp,order.max=order.max) },
										 wts,df,ncp,pow,SIMPLIFY=FALSE)
	kappa <- Reduce('+', subkappa)
	return(kappa)
}
sumchipow.support <- function(wts,df,ncp,pow) {
	minv <- ifelse(min(wts) < 0,-Inf,0)
	maxv <- ifelse(max(wts) > 0,Inf,0)
	retval <- c(minv,maxv)
	return(retval)
}

# dsumchipow, psumchipow, qsumchipow, rsumchipow#FOLDUP
#' @title The sum of (non-central) chi-squares raised to powers distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the weighted sum of non-central
#' chi-squares taken to powers.
#'
#' @details
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares, where \eqn{\nu_i}{v_i}
#' are the degrees of freedom, and \eqn{\delta_i}{delta_i} are the
#' non-centrality parameters.  
#' Let \eqn{w_i} and \eqn{p_i} be given constants. Suppose
#' \deqn{Y = \sum_i w_i X_i^{p_i}.}{Y = sum w_i (X_i)^(p_i).}
#' Then \eqn{Y}{Y} follows a weighted sum of chi-squares to power distribution. 
#'
#' @usage
#'
#' dsumchipow(x, wts, df, ncp=0, pow=1, log = FALSE, order.max=6)
#'
#' psumchipow(q, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qsumchipow(p, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rsumchipow(n, wts, df, ncp=0, pow=1)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param wts the vector of weights. 
#' This is recycled against the \code{df, ncp, pow}, but not against the \code{x,q,p,n}.
#' @param df the vector of degrees of freedom. 
#' This is recycled against the \code{wts, ncp, pow}, but not against the \code{x,q,p,n}.
#' @param ncp the vector of non-centrality parameters. 
#' This is recycled against the \code{wts, df, pow}, but not against the \code{x,q,p,n}.
#' @param ncp the vector of the power parameters. 
#' This is recycled against the \code{wts, df, ncp}, but not against the \code{x,q,p,n}.
#'
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#'
#' @return \code{dsumchipow} gives the density, \code{psumchipow} gives the 
#' distribution function, \code{qsumchipow} gives the quantile function, 
#' and \code{rsumchipow} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dsumchipow psumchipow qsumchipow rsumchipow
#' @seealso 
#' The special cases of the sum chi-square distribution
#' \code{\link{dsumchisq}}, and the
#' sum chi distribution \code{\link{dsumchi}}.
#' @template etc
#' @examples 
#' wts <- c(1,-3,4)
#' df <- c(100,20,10)
#' ncp <- c(5,3,1)
#' pow <- c(1,0.5,1)
#' rvs <- rsumchipow(128, wts, df, ncp, pow)
#' dvs <- dsumchipow(rvs, wts, df, ncp, pow)
#' qvs <- psumchipow(rvs, wts, df, ncp, pow)
#' pvs <- qsumchipow(ppoints(length(rvs)), wts, df, ncp, pow)
#' @rdname dsumchipow
#' @name sumchipow
#' @export 
dsumchipow <- function(x, wts, df, ncp=0, pow=1, log = FALSE, order.max=6) {
	kappa <- sumchipow_cumuls(wts,df,ncp,pow,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,support=sumchipow.support(wts,df,ncp,pow),log=log)
	return(retval)
}
#' @export
psumchipow <- function(q, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumchipow_cumuls(wts,df,ncp,pow,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,support=sumchipow.support(wts,df,ncp,pow),lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qsumchipow <- function(p, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumchipow_cumuls(wts,df,ncp,pow,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,support=sumchipow.support(wts,df,ncp,pow),lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rsumchipow <- function(n, wts, df, ncp=0, pow=1) {
	subX <- mapply(function(w,dd,nn,pp) { w * (rchisq(n,df=dd,ncp=nn) ^ pp) },
										 wts,df,ncp,pow,SIMPLIFY=FALSE)
	X <- Reduce('+', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
