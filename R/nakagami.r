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
# Created: 2015.02.08
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# dnakagami, pnakagami, qnakagami, rnakagami#FOLDUP
#' @title The Nakagami distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the Nakagami distribution.
#'
#' @details
#'
#' Let \eqn{Y \sim \mathrm{Gamma}\left(m, \Omega/m\right)}{Y ~ Gamma(m,Omega/m)}.
#' Then \eqn{\sqrt{Y}}{sqrt(Y)} takes a \emph{Nakagami distribution}
#' with shape \eqn{m}{m} and spread \eqn{\Omega}{Omega}.
#' Equivalently, a Nakagami is a rescaled chi random variable.
#'
#' @usage
#'
#' dnakagami(x, m, Omega=1, log = FALSE)
#'
#' pnakagami(q, m, Omega=1, lower.tail = TRUE, log.p = FALSE)
#'
#' qnakagami(p, m, Omega=1, lower.tail = TRUE, log.p = FALSE)
#'
#' rnakagami(n, m, Omega=1)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param m the shape parameter.
#' @param Omega the spread parameter.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#'  \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#'
#' @references
#'
#' M. Nakagami. "The m-Distribution, a general formula of intensity of rapid fading."
#' Statistical Methods in Radio Wave Propagation (1960): pp 3-36. 
#'
#' @keywords distribution 
#' @return \code{dnakagami} gives the density, \code{pnakagami} gives the 
#' distribution function, \code{qnakagami} gives the quantile function, 
#' and \code{rnakagami} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dnakagami pnakagami qnakagami rnakagami
#' @seealso gamma distribution functions, \code{\link{dgamma}, \link{pgamma}, \link{qgamma}, \link{rgamma}},
#' the chi distribution functions, \code{\link{dchi}, \link{pchi}, \link{qchi}, \link{rchi}}
#' @note
#' This is a thin wrapper on the gamma distribution functions.
#' @template etc
#' @examples 
#' rvs <- rnakagami(128, 20, 1)
#' dvs <- dnakagami(rvs, 20, 1)
#' @rdname dnakagami
#' @name nakagami
#' @export 
dnakagami <- function(x,m,Omega=1,log = FALSE) {
	dens <- dgamma(x^2,shape=m,rate=Omega/m,log=log)
	if (log)
		dens <- dens + log(2*x)
	else
		dens <- dens * (2*x)
	return(dens)
}
#' @export
pnakagami <- function(q,m,Omega=1,lower.tail = TRUE,log.p = FALSE) {
	cdf <- pgamma(q^2,shape=m,rate=Omega/m,lower.tail=lower.tail,log.p=log.p)
	return(cdf)
}
#' @export 
qnakagami <- function(p,m,Omega=1,lower.tail = TRUE,log.p = FALSE) {
	qtile <- qgamma(p,shape=m,rate=Omega/m,lower.tail=lower.tail,log.p=log.p)
	qtile <- sqrt(qtile)
	return(qtile)
}
#' @export 
rnakagami <- function(n,m,Omega=1) {
	X <- sqrt(rgamma(n,shape=m,rate=Omega/m))
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
