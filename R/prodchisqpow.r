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
# Created: 2015.03.07
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the moments of the prodchisqpow
# distribution. 
prodchisqpow_moments <- function(df,ncp=0,pow=1,order.max=3) {
	submu <- mapply(function(dd,nn,pp) 
										 { chipow_moms(df=dd,ncp=nn,pow=pp,order.max=order.max) },
										 df,ncp,pow,SIMPLIFY=FALSE)
	mu <- Reduce('*', submu)
	return(mu)
}
# compute the cumulants of the prodchisqpow
# distribution. 
prodchisqpow_cumuls <- function(df,ncp=0,pow=1,order.max=3) {
	kappa <- moment2cumulant(prodchisqpow_moments(df,ncp,pow,order.max))
	return(kappa)
}

# dprodchisqpow, pprodchisqpow, qprodchisqpow, rprodchisqpow#FOLDUP
#' @title The product of (non-central) chi-squares raised to powers distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the product of non-central
#' chi-squares taken to powers.
#'
#' @details
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares, where \eqn{\nu_i}{v_i}
#' are the degrees of freedom, and \eqn{\delta_i}{delta_i} are the
#' non-centrality parameters.  
#' Let \eqn{p_i} be given constants. Suppose
#' \deqn{Y = \prod_i X_i^{p_i}.}{Y = prod w_i (X_i)^(p_i).}
#' Then \eqn{Y}{Y} follows a product of chi-squares to power distribution. 
#'
#' @usage
#'
#' dprodchisqpow(x, df, ncp=0, pow=1, log = FALSE, order.max=5)
#'
#' pprodchisqpow(q, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=5)
#'
#' qprodchisqpow(p, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=5)
#'
#' rprodchisqpow(n, df, ncp=0, pow=1)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the vector of degrees of freedom. 
#' This is recycled against the \code{ncp, pow}, but not against the \code{x,q,p,n}.
#' @param ncp the vector of non-centrality parameters. 
#' This is recycled against the \code{df, pow}, but not against the \code{x,q,p,n}.
#' @param pow the vector of the power parameters. 
#' This is recycled against the \code{df, ncp}, but not against the \code{x,q,p,n}.
#'
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#'
#' @return \code{dprodchisqpow} gives the density, \code{pprodchisqpow} gives the 
#' distribution function, \code{qprodchisqpow} gives the quantile function, 
#' and \code{rprodchisqpow} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dprodchisqpow pprodchisqpow qprodchisqpow rprodchisqpow
#' @seealso 
#' The upsilon distribution, 
#' \code{\link{dupsilon},\link{pupsilon}},\code{\link{qupsilon},\link{rupsilon}}.
#' The sum of chi-square powers distribution, 
#' \code{\link{dsumchisqpow},\link{psumchisqpow}},\code{\link{qsumchisqpow},\link{rsumchisqpow}}.
#' @template etc
#' @examples 
#' df <- c(100,20,10)
#' ncp <- c(5,3,1)
#' pow <- c(1,0.5,1)
#' rvs <- rprodchisqpow(128, df, ncp, pow)
#' dvs <- dprodchisqpow(rvs, df, ncp, pow)
#' qvs <- pprodchisqpow(rvs, df, ncp, pow)
#' pvs <- qprodchisqpow(ppoints(length(rvs)), df, ncp, pow)
#' @rdname dprodchisqpow
#' @name prodchisqpow
#' @export 
dprodchisqpow <- function(x, df, ncp=0, pow=1, log = FALSE, order.max=5) {
	kappa <- prodchisqpow_cumuls(df,ncp,pow,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export
pprodchisqpow <- function(q, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=5) {
	kappa <- prodchisqpow_cumuls(df,ncp,pow,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qprodchisqpow <- function(p, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=5) {
	kappa <- prodchisqpow_cumuls(df,ncp,pow,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rprodchisqpow <- function(n, df, ncp=0, pow=1) {
	subX <- mapply(function(dd,nn,pp) { (rchisq(n,df=dd,ncp=nn) ^ pp) },
										 df,ncp,pow,SIMPLIFY=FALSE)
	X <- Reduce('*', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
