# Copyright 2014-2016 Steven E. Pav. All Rights Reserved.
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
# Created: 2016.03.06
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the moments of the prodnormal
# distribution. 
prodnormal_moms <- function(mu,sigma,order.max=3) {
	subiota <- mapply(function(mm,ss) { norm_moms(mm,ss,order.max) },
										 mu,sigma,SIMPLIFY=FALSE)
	iota <- Reduce('*', subiota)
	return(iota)
}
prodnormal_cumuls <- function(mu,sigma,order.max=3) {
	kappa <- PDQutils::moment2cumulant(prodnormal_moms(mu,sigma,order.max))
	return(kappa)
}

# dprodnormal, pprodnormal, qprodnormal, rprodnormal#FOLDUP
#' @title The product of normal random variates.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the product of indepdendent
#' normal random variables.
#'
#' @details
#'
#' Let \eqn{Z_i \sim \mathcal{N}\left(\mu_i, \sigma_i^2\right)}{Z_i ~ N(mu_i, sigma_i^2)}
#' be independently distributed normal variates, with means \eqn{\mu_i}{mu_i}
#' and variances \eqn{\sigma_i^2}{sigma_i^2}.
#' Suppose \deqn{Y = \prod_i Z_i.}{Z = prod (Z_i).}
#' Then \eqn{Y}{Y} follows a product of normals distribution. 
#'
#' @usage
#'
#' dprodnormal(x, mu, sigma, log = FALSE, order.max=5)
#'
#' pprodnormal(q, mu, sigma, lower.tail = TRUE, log.p = FALSE, order.max=5)
#'
#' qprodnormal(p, mu, sigma, lower.tail = TRUE, log.p = FALSE, order.max=5)
#'
#' rprodnormal(n, mu, sigma)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param mu the vector of means.
#' This is recycled against the \code{sigma}, but not against the \code{x,q,p,n}.
#' @param sigma the vector of standard deviations.
#' This is recycled against the \code{mu}, but not against the \code{x,q,p,n}.
#'
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @return \code{dprodnormal} gives the density, \code{pprodnormal} gives the 
#' distribution function, \code{qprodnormal} gives the quantile function, 
#' and \code{rprodnormal} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dprodnormal pprodnormal qprodnormal rprodnormal
#' @examples 
#' mu <- c(100,20,10)
#' sigma <- c(10,50,10)
#' rvs <- rprodnormal(128, mu, sigma)
#' dvs <- dprodnormal(rvs, mu, sigma)
#' qvs <- pprodnormal(rvs, mu, sigma)
#' pvs <- qprodnormal(ppoints(length(rvs)), mu, sigma)
#' @rdname dprodnormal
#' @name prodnormal
#' @export 
dprodnormal <- function(x, mu, sigma, log = FALSE, order.max=5) {
	kappa <- prodnormal_cumuls(mu, sigma, order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export
pprodnormal <- function(q, mu, sigma, lower.tail = TRUE, log.p = FALSE, order.max=5) {
	kappa <- prodnormal_cumuls(mu, sigma, order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qprodnormal <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE, order.max=5) {
	kappa <- prodnormal_cumuls(mu, sigma, order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rprodnormal <- function(n, mu, sigma) {
	subX <- mapply(function(mm,ss) { rnorm(n,mean=mm,sd=ss) },
										 mu,sigma,SIMPLIFY=FALSE)
	X <- Reduce('*', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
