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
# Created: 2015.02.27
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the cumulants of the sumchi
# distribution. 
sumchi.cumuls <- function(wts,df,order.max=3) {
	subkappa <- mapply(function(w,dd) 
										 { (w ^ (1:order.max)) * chi.cumuls(df=dd,order.max=order.max) },
										 wts,df,SIMPLIFY=FALSE)
	kappa <- Reduce('+', subkappa)
	return(kappa)
}

# dsumchi, psumchi, qsumchi, rsumchi#FOLDUP
#' @title The sum of (non-central) chi-squares distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the weighted sum of non-central
#' chi-squares.
#'
#' @details
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares. Let \eqn{a_i} be
#' given constants. Suppose
#' \deqn{Y = \sum_i a_i X_i.}{Y = sum a_i X_i.}
#' Then \eqn{Y}{Y} follows a weighted sum of chi-squares distribution. When
#' the weights are all one, and the chi-squares are all central, then 
#' \eqn{Y}{Y} also follows a chi-square distribution.
#'
#' @usage
#'
#' dsumchi(x, wts, df, order.max=6, log = FALSE)
#'
#' psumchi(q, wts, df, order.max=6, lower.tail = TRUE, log.p = FALSE)
#'
#' qsumchi(p, wts, df, order.max=6, lower.tail = TRUE, log.p = FALSE)
#'
#' rsumchi(n, wts, df)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param wts the vector of weights. We do \emph{not} vectorize over
#' \code{wts}, except against \code{df}.
#' @param df the vector of degrees of freedom. We do \emph{not} vectorize over
#' \code{df}, except against \code{wts}.
#' @template distribution
#' @template apx_distribution
#'
#' @return \code{dsumchi} gives the density, \code{psumchi} gives the 
#' distribution function, \code{qsumchi} gives the quantile function, 
#' and \code{rsumchi} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dsumchi psumchi qsumchi rsumchi
#' @seealso sum of chi-square distribution functions,
#' \code{\link{dsumchisq}, \link{psumchisq}, \link{qsumchisq}, \link{rsumchisq}},
#' @template etc
#' @examples 
#' wts <- c(1,-3,4)
#' df <- c(100,20,10)
#' rvs <- rsumchi(128, wts, df)
#' dvs <- dsumchi(rvs, wts, df)
#' qvs <- psumchi(rvs, wts, df)
#' pvs <- qsumchi(ppoints(length(rvs)), wts, df)
#' @rdname dsumchi
#' @name sumchi
#' @export 
dsumchi <- function(x, wts, df, order.max=6, log = FALSE) {
	kappa <- sumchi.cumuls(wts,df,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export
psumchi <- function(q, wts, df, order.max=6, lower.tail = TRUE, log.p = FALSE) {
	kappa <- sumchi.cumuls(wts,df,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qsumchi <- function(p, wts, df, order.max=6, lower.tail = TRUE, log.p = FALSE) {
	kappa <- sumchi.cumuls(wts,df,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa)
	return(retval)
}
#' @export 
rsumchi <- function(n, wts, df) {
	subX <- mapply(function(w,dd) { w * sqrt(rchisq(n,df=dd)) },
										 wts,df,SIMPLIFY=FALSE)
	X <- Reduce('+', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
