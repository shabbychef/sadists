# Copyright 2014-2014 Steven E. Pav. All Rights Reserved.
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
# Created: 2014.02.21
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# dchi, pchi, qchi, rchi#FOLDUP
#' @title The (non-central) chi distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the (non-central) chi distribution.
#'
#' @details
#'
#' Let \eqn{X_i \sim \mathcal{N}\left(\mu_i,1\right)}{X_i ~ N(u_i,1)} 
#' independently for \eqn{1 \le i \le k}{1 <= i <= k}. The random
#' variable 
#' \deqn{T = \sqrt{\sum_i X_i^2}}{T = sqrt(sum_i X_i^2)}
#' follows a \emph{non-central chi distribution} with parameters
#' \eqn{k, \lambda}{k, lambda}, where \eqn{k} is the degrees
#' of freedom, and 
#' \eqn{\lambda = \sum_i \mu_i^2}{lambda = sum_i u_i^2} is
#' the non-centrality parameter. 
#' Equivalently, a non-central chi is the square root of a non-central
#' chi-square random variable.
#'
#' @usage
#'
#' dchi(x, df, ncp = 0, log = FALSE)
#'
#' pchi(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#'
#' qchi(p, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)
#'
#' rchi(n, df, ncp = 0)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the degrees of freedom, \eqn{k}.
#' @param ncp the non-centrality parameter, \eqn{\lambda}{lambda}.
#' When equal to zero, we recover the chi distribution.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#'  \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @inheritParams dchisq
#' @inheritParams pchisq
#' @inheritParams qchisq
#' @inheritParams rchisq
#'
#' @keywords distribution 
#' @return \code{dchi} gives the density, \code{pchi} gives the 
#' distribution function, \code{qchi} gives the quantile function, 
#' and \code{rchi} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases pchi qchi rchi
#' @seealso t distribution functions, \code{\link{dchisq}, \link{pchisq}, \link{qchisq}, \link{rchisq}}
#' @note
#' This is a thin wrapper on the chi-square distribution functions.
#' @template etc
#' @examples 
#' rvs <- rchi(128, 20, 1)
#' dvs <- dchi(rvs, 20, 1)
#' pvs.H0 <- pchi(rvs, 20, 0)
#' pvs.HA <- pchi(rvs, 20, 1)
#' \dontrun{
#' plot(ecdf(pvs.H0))
#' plot(ecdf(pvs.HA))
#' }
#' @export 
dchi <- function(x,df,ncp = 0,log = FALSE) {
	dens <- dchisq(x^2,df=df,ncp=ncp,log=log)
	if (log)
		dens <- dens + log(2*x)
	else
		dens <- dens * (2*x)
	return(dens)
}
#' @export
pchi <- function(q,df,ncp = 0,lower.tail = TRUE,log.p = FALSE) {
	if (hasArg(ncp))
		cdf <- pchisq(q^2,df=df,ncp=ncp,lower.tail=lower.tail,log.p=log.p)
	else
		cdf <- pchisq(q^2,df=df,lower.tail=lower.tail,log.p=log.p)
	return(cdf)
}
#' @export 
qchi <- function(p,df,ncp = 0,lower.tail = TRUE,log.p = FALSE) {
	if (hasArg(ncp))
		qtile <- qchisq(p,df=df,ncp=ncp,lower.tail=lower.tail,log.p=log.p)
	else
		qtile <- qchisq(p,df=df,lower.tail=lower.tail,log.p=log.p)
	if (log.p)
		qtile <- 0.5 * qtile
	else
		qtile <- sqrt(qtile)
	return(qtile)
}
#' @export 
rchi <- function(n,df,ncp = 0) {
	if (hasArg(ncp))
		X <- sqrt(rchisq(n,df=df,ncp=ncp))
	else
		X <- sqrt(rchisq(n,df=df))
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
