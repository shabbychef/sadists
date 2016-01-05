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
# Created: 2015.06.04
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav


# generate generalized gamma variates:
# call gamma, then take to a power
rgengamma <- function(n,a,d,p) {
	# k, theta are shape, scale
	k <- d / p
	theta <- a ^ p
	retv <- rgamma(n, shape=k, scale=theta)
	retv <- retv ^ (1/p)
	retv
}

# moments of the generalized gamma.
gengamma_moms <- function(a,d,p,order.max=3,orders=1:order.max,log=FALSE) {
	dbyp <- d/p
	if (log) {
		retval <- orders * log(a) + lgamrat((orders/p) + dbyp,dbyp)
	} else {
		retval <- (a^(orders)) * gamrat((orders/p) + dbyp,dbyp)
	}
	return(retval)
}

# the generalized gamma.
gengamma_cumuls <- function(a,d,p,order.max=3) {
	retval <- moment2cumulant(gengamma_moms(a,d,p,order.max=order.max,log=FALSE))
	return(retval)
}

# compute the cumulants of the sumgengamma
# distribution. 
sumgengamma_cumuls <- function(wts,a,d,pow,order.max=3) {
	subkappa <- mapply(function(ww,aa,dd,pp) 
										 { (ww ^ (1:order.max)) * gengamma_moms(a=aa,d=dd,p=pp,order.max=order.max) },
										 wts,a,d,pow,SIMPLIFY=FALSE)
	kappa <- Reduce('+', subkappa)
	return(kappa)
}
sumgengamma_support <- function(wts,a,d,pow) {
	minv <- ifelse(min(wts) < 0,-Inf,0)
	maxv <- ifelse(max(wts) > 0,Inf,0)
	retval <- c(minv,maxv)
	return(retval)
}

# dsumgengamma, psumgengamma, qsumgengamma, rsumgengamma#FOLDUP
#' @title The sum of generalized gammas distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the weighted sum of gamma
#' variates.
#'
#' @details
#'
#' Let \eqn{X_i \sim G(a_i, d_i, p_i)}{X_i ~ G(a_i,d_i,p_i)}
#' be independently distributed generalized gamma variates with parameters
#' \eqn{a_i}{a_i}, \eqn{d_i}{d_i}, and inverse power \eqn{p_i}{p_i}.
#' Let \eqn{w_i} be given constants. Suppose
#' \deqn{Y = \sum_i w_i X_i.}{Y = sum w_i X_i.}
#' Then \eqn{Y}{Y} follows a weighted sum of generalized gammas distribution. 
#'
#' @usage
#'
#' dsumgengamma(x, wts, a, d, pow, log = FALSE, order.max=6)
#'
#' psumgengamma(q, wts, a, d, pow, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qsumgengamma(p, wts, a, d, pow, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rsumgengamma(n, wts, a, d, pow)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param wts the vector of weights. 
#' This is recycled against the \code{a, d, pow}, but not against the \code{x,q,p,n}.
#' @param a the vector of the 'a' parameters of the variates.
#' This is recycled against the \code{wts, d, pow}, but not against the \code{x,q,p,n}.
#' @param d the vector of the 'd' parameters of the variates.
#' This is recycled against the \code{wts, a, pow}, but not against the \code{x,q,p,n}.
#' @param pow the vector of the 'p' parameters of the variates. 
#' This is recycled against the \code{wts, a, d}, but not against the \code{x,q,p,n}.
#'
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @template ref-gengamma
#'
#' @return \code{dsumgengamma} gives the density, \code{psumgengamma} gives the 
#' distribution function, \code{qsumgengamma} gives the quantile function, 
#' and \code{rsumgengamma} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dsumgengamma psumgengamma qsumgengamma rsumgengamma
#' @seealso 
#' The sum of chi squares to powers distribution,
#' \code{\link{dsumchisqpow},\link{psumchisqpow},\link{qsumchisqpow},\link{rsumchisqpow}}.
#' @template etc
#' @examples 
#' wts <- c(1,-3,4)
#' a <- c(10,2,5)
#' d <- c(5,30,1)
#' pow <- c(1,2,1)
#' rvs <- rsumgengamma(128, wts, a, d, pow)
#' dvs <- dsumgengamma(rvs, wts, a, d, pow)
#' qvs <- psumgengamma(rvs, wts, a, d, pow)
#' pvs <- qsumgengamma(ppoints(length(rvs)), wts, a, d, pow)
#' @rdname dsumgengamma
#' @name sumgengamma
#' @export 
dsumgengamma <- function(x, wts, a, d, pow, log = FALSE, order.max=6) {
	kappa <- sumgengamma_cumuls(wts,a,d,pow,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,support=sumgengamma_support(wts,a,d,pow),log=log)
	return(retval)
}
#' @export
psumgengamma <- function(q, wts, a, d, pow, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumgengamma_cumuls(wts,a,d,pow,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,support=sumgengamma_support(wts,a,d,pow),
																		 lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qsumgengamma <- function(p, wts, a, d, pow, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumgengamma_cumuls(wts,a,d,pow,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,support=sumgengamma_support(wts,a,d,pow),
															lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rsumgengamma <- function(n,wts,a,d,pow) {
	subX <- mapply(function(ww,aa,dd,pp) { ww * (rgengamma(n,a=aa,d=dd,p=pp)) },
										 wts,a,d,pow,SIMPLIFY=FALSE)
	X <- Reduce('+', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
