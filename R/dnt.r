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

# Created: 2014.02.14
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# ddnt, pdnt, qdnt, rdnt#FOLDUP
#' @title The doubly non-central t distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the doubly non-central t-distribution.
#'
#' @details
#'
#' Let \eqn{X \sim \mathcal{N}\left(\mu,1\right)}{X ~ N(u,1)} independently
#' of \eqn{Y \sim \chi^2\left(k,\theta\right)}{Y ~ x^2(k,theta)}. The 
#' random variable
#' \deqn{T = \frac{X}{\sqrt{Y/k}}}{T = X / sqrt(Y/k)}
#' takes a \emph{doubly non-central t-distribution} with parameters
#' \eqn{k, \mu, \theta}{k, mu, theta}.
#'
#' @usage
#'
#' ddnt(x, k, mu, theta, log = FALSE)
#'
#' pdnt(q, k, mu, theta, lower.tail = TRUE, log.p = FALSE)
#'
#' qdnt(p, k, mu, theta, lower.tail = TRUE, log.p = FALSE)
#'
#' rdnt(n, k, mu, theta)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param k the degrees of freedom in the denominator.
#' @param mu the numerator non-centrality parameter, \eqn{\mu}{mu}.
#' @param theta the denominator non-centrality parameter, \eqn{\theta}{theta}.
#' When equal to zero, we recover the singly non-central t-distribution.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#'  \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @inheritParams dt
#' @inheritParams pt
#' @inheritParams qt
#' @inheritParams rt
#'
#' @keywords distribution 
#' @return \code{ddnt} gives the density, \code{pdnt} gives the 
#' distribution function, \code{qdnt} gives the quantile function, 
#' and \code{rdnt} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases pdnt qdnt rdnt
#' @seealso t-distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
#' @note
#' This is a thin wrapper on the t distribution. 
#' The functions \code{\link{dt}, \link{pt}, \link{qt}} can accept ncp from
#' limited range (\eqn{|\delta|\le 37.62}{delta <= 37.62}). Some corrections
#' may have to be made here for large \code{zeta}.
#' @template etc
#' @template dnt
#' @examples 
#' rvs <- rdnt(128, 20, 1, 1)
#' dvs <- ddnt(rvs, 20, 1, 1)
#' pvs.H0 <- pdnt(rvs, 20, 0, 1)
#' pvs.HA <- pdnt(rvs, 20, 1, 1)
#' \dontrun{
#' plot(ecdf(pvs.H0))
#' plot(ecdf(pvs.HA))
#' }
#' # compare to singly non-central
#' dv1 <- ddnt(1, k=10, mu=5, theta=0, log=FALSE)
#' dv2 <- dt(1, df=10, ncp=5, log=FALSE)
#' pv1 <- pdnt(1, k=10, mu=5, theta=0, log.p=FALSE)
#' pv11 <- pdnt(1, k=10, mu=5, theta=0.001, log.p=FALSE)
#' v2 <- pt(1, df=10, ncp=5, log.p=FALSE)
#'
#' q1 <- qdnt(pv1, k=10, mu=5, theta=0, log.p=FALSE)
#'
# listing 10.13
.ddnt <- function(x, k, mu, theta, log=FALSE) {
	aterm <- x * mu * sqrt(2/k)
	kon <- (-(theta + mu^2)/2) - log(pi*k)/2;
	w <- (1+(x^2/k))
	logw <- log(w)

	doubnoncentterm <- function(j) {
		if (aterm == 0) {
			if (j == 0)
				term1 <- 0
			else 
				return(0)
		} else 
			term1 <- j * log(aterm)
		z <- (k+j+1)/2
		term <- kon + term1 + lgamma(z) - lgamma(j+1) - lgamma(k/2) - z * logw
		f <- Re(hypergeo::genhypergeo(z,k/2,theta/(2*w)))
		y <- Re(exp(term + log(f)))
	}
	kum <- 0
	for (jjj in c(0:300)) {
		kum <- kum + doubnoncentterm(jjj)
	}
	# sorry, no better way to do this.
	if (log)
		kum <- log(kum)
	return(kum)
}
#' @export 
ddnt <- Vectorize(.ddnt,
									vectorize.args = c("x","k","mu","theta"),
									SIMPLIFY = TRUE)
#' listing 10.15
.pdnt <- function(q, k, mu, theta, lower.tail = TRUE, log.p = FALSE) {
	if (theta == 0)
		return(pt(q,df=k,ncp=mu,lower.tail=lower.tail,log.p=log.p))
	partsum <- 0
	inc <- 50
	lo <- 0
	hi <- inc-1
	done <- FALSE
	while (!done) {
		ivec <- lo:hi
		lnw <- (-theta/2) + ivec * log(theta/2) - lgamma(ivec+1)
		k2i <- k + 2*ivec
		logcdfpart <- pt(q * sqrt(k2i/k),df = k2i,ncp=mu,log.p=TRUE)
		newpartvec <- exp(lnw + logcdfpart)
		newpartsum <- sum(newpartvec)
		partsum <- partsum + newpartsum
		lo <- hi+1
		hi <- hi+inc
		done <- (newpartvec[1] >= newpartvec[length(newpartvec)]) && (newpartsum < 1e-12)
	}
	cdf <- Re(partsum)
	if (! lower.tail)
		cdf <- 1 - cdf
	if (log.p)
		cdf <- log(cdf)
	return(cdf)
}
#' @export 
pdnt <- Vectorize(.pdnt,
									vectorize.args = c("q","k","mu","theta"),
									SIMPLIFY = TRUE)
#' uh, invert it? numerically?
.qdnt <- function(p, k, mu, theta, lower.tail = TRUE, log.p = FALSE) {
	if (lower.tail) 
		zerf <- function(q) {
			.pdnt(q,k,mu,theta,lower.tail=lower.tail,log.p=log.p) - p
		}
	else
		zerf <- function(q) {
			p - .pdnt(q,k,mu,theta,lower.tail=lower.tail,log.p=log.p)
		}
	q0 <- qt(p, df=k, ncp=mu, lower.tail=lower.tail, log.p=log.p)
	v0 <- zerf(q0)
	if (v0 == 0) {
		return(q0)
	} else if (v0 > 0) {
		ub <- q0
		fub <- v0
		lb <- q0 - 0.1 * max(0.01,abs(q0))
		flb <- zerf(lb)
# 2FIX: beware infs!
		while (flb > 0) {
			# drag ub down too
			ub <- lb
			fub <- flb
			lb <- q0 - 2 * (q0 - lb)
			flb <- zerf(lb)
		}
	} else {
		lb <- q0
		flb <- v0
		ub <- q0 + 0.1 * max(0.01,abs(q0))
		fub <- zerf(ub)
# 2FIX: beware infs!
		while (fub < 0) {
			# drag lb up too
			lb <- ub
			flb <- fub
			ub <- q0 + 2 * (ub - q0)
			fub <- zerf(ub)
		}
	}
	rfnd <- uniroot(zerf,c(lb,ub),f.lower=flb,f.upper=fub)
	return(rfnd$root)
}
#' @export 
qdnt <- Vectorize(.qdnt,
									vectorize.args = c("p","k","mu","theta"),
									SIMPLIFY = TRUE)
#' @export 
rdnt <- function(n, k, mu, theta) {
	X <- rnorm(n,mean=mu,sd=1)
	Y <- rchisq(n,df=k,ncp=theta)
	Z <- X / sqrt(Y / k)
	return(Z)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
