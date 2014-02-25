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

# dkprime, pkprime, qkprime, rkprime#FOLDUP
#' @title The K prime distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the K prime distribution.
#'
#' @details
#'
#' Suppose \eqn{y \sim \chi^2\left(\nu_1\right)}{y ~ x^2(v1)}, and
#' \eqn{x \sim t \left(\nu_2, a\sqrt{y/\nu_1}/b\right)}{x ~ t(v2,(a/b) sqrt(y/v1))}.
#' Then the random variable
#' \deqn{T = b x}{T = b x}
#' takes a K prime distribution with parameters 
#' \eqn{\nu_1, \nu_2, a, b}{v1, v2, a, b}. In Lecoutre's terminology,
#' \eqn{T \sim K'_{\nu_1, \nu_2}\left(a, b\right)}{T ~ K'_v1,v2(a,b)}
#'
#' Equivalently, we can think of
#' \deqn{T = \frac{b Z + a \sqrt{\chi^2_{\nu_1} / \nu_1}}{\sqrt{\chi^2_{\nu_2} / \nu_2}}}{T = (bZ + a sqrt(chi2_v1/v1)) / sqrt(chi2_v2/v2)}
#' where \eqn{Z} is a standard normal, and the normal and the (central) chi-squares are
#' independent of each other. When \eqn{a=0}{a=0} we recover
#' a central t distribution; 
#' when \eqn{\nu_1=\infty}{v1=inf} we recover a rescaled non-central t distribution;
#' when \eqn{b=0}{b=0}, we get a rescaled square root of a central F
#' distribution; when \eqn{\nu_2=\infty}{v2=inf}, we recover a 
#' Lambda prime distribution.
#'
#' @usage
#'
#' dkprime(x, v1, v2, a, b = 1, log = FALSE)
#'
#' pkprime(q, v1, v2, a, b = 1, lower.tail = TRUE, log.p = FALSE)
#'
#' qkprime(p, v1, v2, a, b = 1, lower.tail = TRUE, log.p = FALSE)
#'
#' rkprime(n, v1, v2, a, b = 1)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param v1 the degrees of freedom in the numerator chisquare. When
#' (positive) infinite, we recover a non-central t 
#' distribution with \code{v2} degrees of freedom and non-centrality
#' parameter \code{a}, scaled by \code{b}.
#' @param v2 the degrees of freedom in the denominator chisquare.
#' When equal to infinity, we recover the Lambda prime distribution.
#' @param a the non-centrality scaling parameter. When equal to zero,
#' we recover the (central) t distribution.
#' @param b the scaling parameter.
#'
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
#' @return \code{dkprime} gives the density, \code{pkprime} gives the 
#' distribution function, \code{qkprime} gives the quantile function, 
#' and \code{rkprime} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dkprime pkprime qkprime rkprime
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}},
#' K square distribution functions, \code{\link{dksquare}, \link{pksquare}, \link{qksquare}, \link{rksquare}},
#' lambda prime distribution functions, \code{\link{dlambdap}, \link{plambdap}, \link{qlambdap}, \link{rlambdap}},
#' @template etc
#' @template ref-kprime
#' @examples 
#' d1 <- dkprime(1, 50, 20, a=0.01)
#' d2 <- dkprime(1, 50, 20, a=0.0001)
#' d3 <- dkprime(1, 50, 20, a=0)
#' d4 <- dkprime(1, 10000, 20, a=1)
#' d5 <- dkprime(1, Inf, 20, a=1)
#'
#' \dontrun{
#' avals <- 10^(seq(1,7,length.out=101))
# 'dvals <- dkprime(1, v1=avals, v2=20, a=1)
#' plot(log10(avals),dvals) 
#' }
#' @rdname dkprime
#' @name kprime
.dkprime <- function(x, v1, v2, a, b = 1, log = FALSE) {
#2FIX: check sane values of v1, v2, a, b?

	# first scale out b;
	x <- x / b
	if (is.infinite(v1)) {
		dens <- dt(x, df=v2, ncp=a, log = log)
	} else if (a == 0) {
		dens <- dt(x, df=v2, log = log)
	} else {
#2FIX: what if is.infinite(v2) ???
		polyterm <- (a*x) * sqrt(4 / ((v1+a^2) * (v2+x^2)))

		f.accum <- function(idx) {
			retval <- ((polyterm/2)^2) * (v2+idx) * (v1+1+idx) / ((idx+2) * (idx+1))
			retval <- cumprod(retval)
		}

		ldenom.0 <- lgamma(v1/2) + lgamma((1+v2)/2)
		ldenom.1 <- lgamma((1+v1)/2) + lgamma((2+v2)/2)

		a0 <- 1
		a1 <- polyterm
		a.even <- a0
		a.odd <- a1
		sum.even <- a0
		sum.odd <- a1

		idx.even <- c(0)
		idx.odd <- c(1)

		ntak <- 100
		done <- FALSE
		while (!done) {
			idx.even <- seq(2 + idx.even[length(idx.even)],by=2,length.out=ntak)
			idx.odd <- seq(2 + idx.odd[length(idx.odd)],by=2,length.out=ntak)
			a.even <- a.even[length(a.even)] * f.accum(idx.even)
			a.odd <- a.odd[length(a.odd)] * f.accum(idx.odd)
			sum.even <- sum.even + sum(a.even)
			sum.odd <- sum.odd + sum(a.odd)
			done <- (idx.even[1] > 100000)
			done <- done || 
				((a.even[length(a.even)] < a.even[1]) && (a.odd[length(a.odd)] < a.odd[1]))
		}

		ldrag <- lgamma(v1/2) + lgamma(v2/2)
		ldenom.0 <- ldenom.0 - ldrag
		ldenom.1 <- ldenom.1 - ldrag

		if (log) {
			cterm <- -0.5 * log(pi * v2)
			cterm <- cterm + (v1/2) * (log(v2) - log(v2 + a^2))
			cterm <- cterm + ((1 + v2)/2) * (log(v2) - log(v2 + x^2))

			# crappy. want a better log expansion..
			proto.dens <- exp(ldenom.0) * sum.even + exp(ldenom.1) * sum.odd

			dens <- cterm + log(proto.dens)
		} else {
			#cterm <- 1 / (sqrt(pi * v2) * exp(ldrag))
			cterm <- 1 / (sqrt(pi * v2))
			cterm <- cterm * (v1/(v2 + a^2))^(v2/2)
			cterm <- cterm * (v2/(v2 + x^2))^((1+v2)/2)

			# crappy. want a better log expansion..
			proto.dens <- exp(ldenom.0) * sum.even + exp(ldenom.1) * sum.odd

			dens <- cterm * proto.dens
		}
	}
	# don't forget the b, though:
	if (log) 
		dens <- dens - log(b)
	else
		dens <- dens / b
	return(dens)
}
#' @export 
dkprime <- Vectorize(.dkprime,
									vectorize.args = c("x","v1","v2","a","b"),
									SIMPLIFY = TRUE)
.pkprime <- function(q, v1, v2, a, b = 1, lower.tail = TRUE, log.p = FALSE) {
}

#' @export 
pkprime <- Vectorize(.pkprime,
									vectorize.args = c("q","v1","v2","a","b"),
									SIMPLIFY = TRUE)
# uh, invert it? numerically?
.qkprime <- function(p, v1, v2, a, b = 1, lower.tail = TRUE, log.p = FALSE) {
}
#' @export 
qkprime <- Vectorize(.qkprime,
									vectorize.args = c("p","v1","v2","a","b"),
									SIMPLIFY = TRUE)
#' @export 
rkprime <- function(n, v1, v2, a, b = 1) {
#2FIX: check for b = 0...
	y <- rchisq(n,df=v1) 
	ncp <- sqrt(y/v1) * (a/b)
	X <- b * rt(n,df=v2,ncp=ncp)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
