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

# Created: 2015.02.08
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <>
# Comments: Steven E. Pav

# compute the 1 through order.max raw, uncentered moment
# of the normal distribution with given mean and standard
# deviation
norm.moms <- function(mu=0,sigma=1,order.max=3) {
	retval <- rep(1,order.max)
	hermi <- hermite.he.polynomials(order.max, normalized=FALSE)
	for (iii in c(1:order.max)) {
		cvals <- abs(coefficients(hermi[[iii+1]]))
		lvals <- mu^(0:iii) * sigma^(iii - (0:iii))
		retval[iii] <- sum(cvals * lvals)
	}
	return(retval)
}

# compute the 1 through order.max raw, uncentered moment
# of the (central) chi distribution with df d.f.
chi.moms <- function(df,order.max=3) {
	jvals <- 1:order.max
	retval <- (2^(jvals/2)) * sapply(jvals,function(j) { gamrat((df+j)/2,df/2) })
	return(retval)
}
# something like a nakagami, but really a scaled chi
schi.moms <- function(df,scal=1,order.max=3) {
	retval <- chi.moms(df=df,order.max=order.max)
	jvals <- 1:order.max
	retval <- retval * ((scal / sqrt(df))^jvals)
	return(retval)
}

# compute the cumulants of the multiple lambda
# prime distribution. this is distributed as
#
# sum_i t_i sqrt(chi^2(df_i) / df_i) + Z
#
# where the chi^2 are independent chi-square
# independent of Z
mlp.cumulants <- function(df,t,order.max=3) {
	# first the Z
	Zmom <- norm.moms(0,1,order.max)
	retval <- moment2cumulant(Zmom)
	for (iii in (1:length(coef))) {
		nmom <- schi.moms(df=df[iii],scal=t[iii],order.max=order.max)
		retval <- retval + moment2cumulant(nmom)
	}
	return(retval)
}

# dmlp, pmlp, qmlp, rmlp#FOLDUP
#' @title The multiple lambda prime distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the 'multiple' lambda prime distribution.
#'
#' @details
#'
#' Suppose \eqn{y_i \sim \chi^2\left(\nu_i\right)}{y_i ~ x^2(v_i)}
#' independently and independently of \eqn{Z}{Z}, a standard normal. 
#' Then 
#' \deqn{T = Z + \sum_i t_i \sqrt{y_i/\nu_i}}{T = Z + sum_i t_i sqrt(y_i/v_i)}
#' takes a 'multiple' lambda prime distribution with parameter vectors
#' \eqn{<\nu_1, \nu_2, \ldots, \nu_k>, <t_1, t_2, ..., t_k>}{<v_1, v_2, ..., v_k>, <t_1, t_2, ..., t_k>}.
#'
#' The multiple lambda prime distribution is used in certain tests of
#' the Sharpe ratio for independent observations, and generalizes
#' the lambda prime random variable.
#'
#' @usage
#'
#' dmlp(x, df, t, log = FALSE)
#'
#' pmlp(q, df, t, lower.tail = TRUE, log.p = FALSE)
#'
#' qmlp(p, df, t, lower.tail = TRUE, log.p = FALSE)
#'
#' rmlp(n, df, t)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the degrees of freedom in the chi square. a vector. we do
#' \emph{not} vectorize over this variable.
#' @param t the scaling parameter on the chi. a vector. should be the same
#' length as \code{df}. we do \emph{not} vectorize over this variable.
#'
#' @return \code{dmlp} gives the density, \code{pmlp} gives the 
#' distribution function, \code{qmlp} gives the quantile function, 
#' and \code{rmlp} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dmlp pmlp qmlp rmlp
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}},
#' K prime distribution functions, \code{\link{dkprime}, \link{pkprime}, \link{qkprime}, \link{rkprime}},
#' @template etc
#' @template distribution
#' @template ref-lambdap
#' @examples 
#' d1 <- dmlp(1, 50, t=0.01)
#' @name mlp
#' @rdname dmlp
#' @export 
dmlp <- function(x, df, t, log = FALSE, order.max=10) {
	kappa <- mlp.cumulants(df,t,order.max=order.max)
	mu.raw <- cumulant2moment(kappa)
	retval <- dapx.gca(x,mu.raw,log=log)
	return(retval)
}
#' @export 
pmlp <- function(q, df, t, lower.tail = TRUE, log.p = FALSE, order.max=10) {
	kappa <- mlp.cumulants(df,t,order.max=order.max)
	mu.raw <- cumulant2moment(kappa)
	retval <- papx.gca(q,mu.raw,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qmlp <- function(p, df, t, lower.tail = TRUE, log.p = FALSE, order.max=10) {
	kappa <- mlp.cumulants(df,t,order.max=order.max)
	retval <- qapx.cf(p,kappa)
	return(retval)
}
#' @export 
rmlp <- function(n, df, t) {
## 2FIX: implement this
	#y <- rchisq(n,df=df) 
	#X <- rnorm(n,mean=t * sqrt(y/df))
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
