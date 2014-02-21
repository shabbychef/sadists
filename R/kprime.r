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
#' Suppose \eqn{y \sim \chi^2\left(\nu_2\right)}{y ~ x^2(v2)}, and
#' \eqn{x \sim t \left(\nu_1, a\sqrt{y/\nu_2}/b\right){x ~ t(v1,(a/b) sqrt(y/v2))}.
#' Then the random variable
#' \deqn{Z = b x}{Z = b x}
#' takes a K prime distribution with parameters 
#' \eqn{\nu_2, \nu_1, a, b}{v2, v1, a, b}. In Lecoutre's terminology,
#' \eqn{Z \sim K'_{\nu_2, \nu_1}\left(a, b\right)}{Z ~ K'_v2,v1(a,b)}
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
#' @param v1 the degrees of freedom in the chisquare which makes up the
#' non-centrality parameter. When equal to infinity, we recover
#' the Lambda prime distribution.
#' @param v2 the degrees of freedom in the non-central t.
#' When equal to infinity, we recover the (scaled) non-central t distribution
#' with \code{v1} degrees of freedom and non-centrality parameter \code{a},
#' scaled by \code{b}.
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
#' @aliases pkprime qkprime rkprime
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
#' K square distribution functions, \code{\link{dksquare}, \link{pksquare}, \link{qksquare}, \link{rksquare}}
#' @template etc
#' @template kprime
#' @examples 
#' d1 <- dkprime(1, 20, 50, a=0.01)
#' d2 <- dkprime(1, 20, 50, a=0.0001)
#' d3 <- dkprime(1, 20, 50, a=0)
#' d4 <- dkprime(1, 20, 10000, a=1)
#' d5 <- dkprime(1, 20, Inf, a=1)
.dkprime <- function(x, v1, v2, a, b = 1, log = FALSE) {
#2FIX: check sane values of v1, v2, a, b?

	# first scale out b;
	x <- x / b
	if (is.infinite(v2)) {
		dens <- dt(x, df=v1, ncp=a, log = log)
	} else if (a == 0) {
		dens <- dt(x, df=v1, log = log)
	} else {
#2FIX: what if is.infinite(v1) ???
		polyterm <- (a*x) * sqrt(4 / ((v2+a^2) * (v1+x^2)))

		f.accum <- function(idx) {
			retval <- ((polyterm/2)^2) * (v1+idx) * (v2+1+idx) / ((idx+2) * (idx+1))
			retval <- cumprod(retval)
		}

		ldenom.0 <- lgamma(v2/2) + lgamma((1+v1)/2)
		ldenom.1 <- lgamma((1+v2)/2) + lgamma((2+v1)/2)

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

		ldrag <- lgamma(v2/2) + lgamma(v1/2)
		ldenom.0 <- ldenom.0 - ldrag
		ldenom.1 <- ldenom.1 - ldrag

		if (log) {
			cterm <- -0.5 * log(pi * v1)
			cterm <- cterm + (v2/2) * (log(v2) - log(v2 + a^2))
			cterm <- cterm + ((1 + v1)/2) * (log(v1) - log(v1 + x^2))

			# crappy. want a better log expansion..
			proto.dens <- exp(ldenom.0) * sum.even + exp(ldenom.1) * sum.odd

			dens <- cterm + log(proto.dens)
		} else {
			#cterm <- 1 / (sqrt(pi * v1) * exp(ldrag))
			cterm <- 1 / (sqrt(pi * v1))
			cterm <- cterm * (v2/(v2 + a^2))^(v2/2)
			cterm <- cterm * (v1/(v1 + x^2))^((1+v1)/2)

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
#'d1 <- dkprime(1, 20, 50, a=0.01)
#'d2 <- dkprime(1, 20, 50, a=0.0001)
#'d3 <- dkprime(1, 20, 50, a=0)
#'d4 <- dkprime(1, 20, 10000, a=1)
#'d4 <- dkprime(1, 20, 1000000, a=1)
#'d5 <- dkprime(1, 20, Inf, a=1)
#'
#'avals <- 10^(seq(1,7,length.out=101))
#'dvals <- dkprime(1, 20, v2=avals, a=1)
#'plot(log10(avals),dvals) 

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
	y <- rchisq(n,df=v2) 
	ncp <- sqrt(y/v2) * (a/b)
	X <- b * rt(n,df=v1,ncp=ncp)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
