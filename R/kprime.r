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
.dkprime <- function(x, v1, v2, a, b = 1, log = FALSE) {
	# first scale out b;
	x <- x / b
	if (is.infinite(v2)) {
		dens <- dt(x, df=v1, ncp=a, log = log)
	} else {
		if (log) {
			cterm <- -0.5 * log(pi * v1) - lgamma(v2 / 2) - lgamma(v1 / 2)
			cterm <- cterm + (v2/2) * (log(v2) - log(v2 + a^2))
			cterm <- cterm + ((1 + v1)/2) * (log(v1) - log(v1 + x^2))
			polyterm <- (a*x) * sqrt(4 / ((v2+a^2) * (v1+x^2)))

		} else {
			cterm <- (1 / (sqrt(pi * v1) * gamma(v2/2) * gamma(v1/2)))
			cterm <- cterm * (v2/(v2 + a^2))^(v2/2)
			cterm <- cterm * (v1/(v1 + x^2))^((1+v1)/2)

			polyterm <- (a*x) * sqrt(4 / ((v2+a^2) * (v1+x^2)))
			a0 <- gamma(v2/2) * gamma((1+v1)/2) 
			a1 <- gamma((1+v2)/2) * gamma((2+v1)/2) * polyterm

			summ <- a0 + a1
			...

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
	y <- rchisq(n,df=v2) 
	ncp <- sqrt(y/v2) * (a/b)
	X <- b * rt(n,df=v1,ncp=ncp)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
