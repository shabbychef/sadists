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

# Created: 2014.02.21
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# dlambdap, plambdap, qlambdap, rlambdap#FOLDUP
#' @title The lambda prime distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the lambda prime distribution.
#'
#' @details
#'
#' Suppose \eqn{y \sim \chi^2\left(\nu\right)}{y ~ x^2(v)}, and
#' \eqn{Z}{Z} is a standard normal. 
#' \deqn{T = Z + t \sqrt{y/\nu}}{T = Z + t sqrt(y/v)}
#' takes a lambda prime distribution with parameters 
#' \eqn{\nu, t}{v, t}.
#' A lambda prime random variable can be viewed as a confidence
#' variable on a non-central t because 
#' \deqn{t = \frac{Z' + T}{\sqrt{y/\nu}}}{t = (Z' + T)/sqrt(y/v)}
#'
#' @usage
#'
#' dlambdap(x, df, t, log = FALSE)
#'
#' plambdap(q, df, t, lower.tail = TRUE, log.p = FALSE)
#'
#' qlambdap(p, df, t, lower.tail = TRUE, log.p = FALSE)
#'
#' rlambdap(n, df, t)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the degrees of freedom in the chi square. 
#' @param t the scaling parameter on the chi.
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
#' @return \code{dlambdap} gives the density, \code{plambdap} gives the 
#' distribution function, \code{qlambdap} gives the quantile function, 
#' and \code{rlambdap} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases plambdap qlambdap rlambdap
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}},
#' K prime distribution functions, \code{\link{dkprime}, \link{pkprime}, \link{qkprime}, \link{rkprime}},
#' @template etc
#' @template ref-lambdap
#' @examples 
#' d1 <- dlambdap(1, 50, t=0.01)
.dlambdap <- function(x, df, t, log = FALSE) {

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
dlambdap <- Vectorize(.dlambdap,
									vectorize.args = c("x","df","t"),
									SIMPLIFY = TRUE)
.plambdap <- function(q, df, t, lower.tail = TRUE, log.p = FALSE) {
}

#' @export 
plambdap <- Vectorize(.plambdap,
									vectorize.args = c("q","df","t"),
									SIMPLIFY = TRUE)
# uh, invert it? numerically?
.qlambdap <- function(p, df, t, lower.tail = TRUE, log.p = FALSE) {
}

#' @export 
qlambdap <- Vectorize(.qlambdap,
									vectorize.args = c("p","df","t"),
									SIMPLIFY = TRUE)
#' @export 
rlambdap <- function(n, df, t) {
	y <- rchisq(n,df=df) 
	X <- rnorm(n,mean=t * sqrt(y/df))
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
