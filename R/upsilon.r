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
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

require(orthopolynom)


# compute the cumulants of the upsilon
# distribution. this is distributed as
#
# sum_i t_i sqrt(chi^2(df_i) / df_i) + Z
#
# where the chi^2 are independent chi-square
# independent of Z
upsilon.cumuls <- function(df,t,order.max=3) {
	# first the Z
	retval <- norm.cumuls(0,1,order.max)
	for (iii in (1:length(coef))) {
		retval <- retval + schi.cumuls(df=df[iii],scal=t[iii],order.max=order.max)
	}
	return(retval)
}

# dupsilon, pupsilon, qupsilon, rupsilon#FOLDUP
#' @title The upsilon distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the upsilon distribution.
#'
#' @details
#'
#' Suppose \eqn{x_i \sim \chi^2\left(\nu_i\right)}{x_i ~ X^2(v_i)}
#' independently and independently of \eqn{Z}{Z}, a standard normal. 
#' Then 
#' \deqn{\Upsilon = Z + \sum_i t_i \sqrt{x_i/\nu_i}}{Y = Z + sum_i t_i sqrt(x_i/v_i)}
#' takes an upsilon distribution with parameter vectors
#' \eqn{[\nu_1, \nu_2, \ldots, \nu_k]', [t_1, t_2, ..., t_k]'}{<v_1, v_2, ..., v_k>, <t_1, t_2, ..., t_k>}.
#'
#' The upsilon distribution is used in certain tests of
#' the Sharpe ratio for independent observations, and generalizes
#' the lambda prime distribution, which can be written as
#' \eqn{Z + t \sqrt{x/\nu}}{Z + t sqrt(x/v)}.
#'
#' @usage
#'
#' dupsilon(x, df, t, log = FALSE, order.max=5)
#'
#' pupsilon(q, df, t, lower.tail = TRUE, log.p = FALSE, order.max=5)
#'
#' qupsilon(p, df, t, lower.tail = TRUE, log.p = FALSE, order.max=5)
#'
#' rupsilon(n, df, t)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the degrees of freedom in the chi square. a vector. we do
#' \emph{not} vectorize over this variable.
#' @param t the scaling parameter on the chi. a vector. should be the same
#' length as \code{df}. we do \emph{not} vectorize over this variable.
#'
#' @return \code{dupsilon} gives the density, \code{pupsilon} gives the 
#' distribution function, \code{qupsilon} gives the quantile function, 
#' and \code{rupsilon} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @note the PDF and CDF are approximated by an Edgeworth expansion; the
#' quantile function is approximated by a Cornish-Fisher expansion.
#' @aliases dupsilon pupsilon qupsilon rupsilon
#' @seealso lambda-prime distribution functions, \code{\link{dlambdap}, \link{plambdap}, \link{qlambdap}, \link{rlambdap}}.
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template ref-lambdap
#' @examples 
#' mydf <- c(100,30,50)
#' myt <- c(-1,3,5)
#' rv <- rupsilon(500, df=mydf, t=myt)
#' d1 <- dupsilon(rv, df=mydf, t=myt)
#' \dontrun{
#' plot(rv,d1)
#' }
#' p1 <- pupsilon(rv, df=mydf, t=myt)
#' # should be nearly uniform:
#' \dontrun{
#' plot(ecdf(p1))
#' }
#' q1 <- qupsilon(ppoints(length(rv)),df=mydf,t=myt)
#' \dontrun{
#' qqplot(x=rv,y=q1)
#' }
#' @name upsilon
#' @rdname dupsilon
#' @export 
dupsilon <- function(x, df, t, log = FALSE, order.max=5) {
	kappa <- upsilon.cumuls(df,t,order.max=order.max)
	#mu.raw <- PDQutils::cumulant2moment(kappa)
	#retval <- PDQutils::dapx_gca(x,mu.raw,log=log)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export 
pupsilon <- function(q, df, t, lower.tail = TRUE, log.p = FALSE, order.max=5) {
	kappa <- upsilon.cumuls(df,t,order.max=order.max)
	#mu.raw <- PDQutils::cumulant2moment(kappa)
	#retval <- PDQutils::papx_gca(q,mu.raw,lower.tail=lower.tail,log.p=log.p)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qupsilon <- function(p, df, t, lower.tail = TRUE, log.p = FALSE, order.max=5) {
	kappa <- upsilon.cumuls(df,t,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa)
	return(retval)
}
#' @export 
rupsilon <- function(n, df, t) {
	retval <- rnorm(n)
	for (iii in (1:length(df))) {
		# or use nakagami?
		retval <- retval + t[iii] * sqrt(rchisq(n,df=df[iii])/df[iii])
	}
	return(retval)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
