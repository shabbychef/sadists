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
#' @inheritParams dt
#' @inheritParams pt
#' @inheritParams qt
#' @inheritParams rt
#'
#' @return \code{dlambdap} gives the density, \code{plambdap} gives the 
#' distribution function, \code{qlambdap} gives the quantile function, 
#' and \code{rlambdap} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dlambdap plambdap qlambdap rlambdap
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}},
#' K prime distribution functions, \code{\link{dkprime}, \link{pkprime}, \link{qkprime}, \link{rkprime}},
#' @template etc
#' @template distribution
#' @template ref-lambdap
#' @examples 
#' d1 <- dlambdap(1, 50, t=0.01)
#' @name lambdap
#' @rdname dlambdap
.dlambdap <- function(x, df, t, log = FALSE) {

	# this _was_ cut and paste from the kprime code,
	# but needs to be written ... 
}
#' @export 
dlambdap <- Vectorize(.dlambdap,
									vectorize.args = c("x","df","t"),
									SIMPLIFY = TRUE)
.plambdap <- function(q, df, t, lower.tail = TRUE, log.p = FALSE) {
# 2FIX: write this.
}

#' @export 
plambdap <- Vectorize(.plambdap,
									vectorize.args = c("q","df","t"),
									SIMPLIFY = TRUE)
# uh, invert it? numerically?
.qlambdap <- function(p, df, t, lower.tail = TRUE, log.p = FALSE) {
# 2FIX: write this.
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
