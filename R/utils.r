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

# Created: 2014.02.23
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the log (gamma(a) / gamma(b))
lgamrat <- function(a,b) {
	return(lgamma(a) - lgamma(b))
}
# compute the ratio gamma(a) / gamma(b)
gamrat <- function(a,b) {
	return(exp(lgamrat(a,b)))
}

# quantile function helper; 
# find the root of a non-decreasing function fnc
# given a single starting point x0.
# starts working outwards from x0 until a sign change is
# detected, then calls uniroot.
#
# the parameters tol and maxiter are passed to uniroot
#
# I think this function is orphaned ... 
uniroot_helper <- function(fnc, x0=0, f0=NULL, xmin=-Inf, xmax=Inf,
													 tol=.Machine$double.eps^0.25, maxiter=1000, 
													 ...) {  # nocov start
	if (is.null(f0)) 
		f0 = fnc(x0,...)
	if (f0 == 0)
		return(x0)
	if (f0 > 0) {
		ub <- x0
		fub <- f0
		lb <- max(xmin,x0 - 0.1 * max(0.01,abs(x0)))
		flb <- fnc(lb,...)
	} else if (f0 < 0) {
		lb <- x0
		flb <- f0
		ub <- min(xmax,x0 + 0.1 * max(0.1,abs(x0)))
		fub <- fnc(ub,...)
	} else
		stop("fnc returned NA or NaN")

	sprod <- sign(flb) * sign(fub)

	while ((sign(flb) * sign(fub) > 0) && (lb > xmin) && (ub < xmax)) {
		if (sign(flb) > 0) {
			# drag ub down too
			ub <- lb
			fub <- flb
			# move lb down
			lb <- max(xmin,x0 - 2 * (x0 - lb))
			flb <- fnc(lb,...)
		} else {
			# drag lb up too
			lb <- ub
			flb <- fub
			ub <- x0 + 2 * (ub - x0)
			fub <- fnc(ub,...)
		}
	}
	# check for xmin and xmax ... 
	if (sign(flb) * sign(fub) > 0) {
		if ((!lb > xmin) && (is.infinite(xmin))) {
			return(xmin)
		} else {
			return(xmax)
		}
	}

	soln <- uniroot(fnc,c(lb,ub),
									f.lower=flb,f.upper=fub,
									tol=tol,maxiter=maxiter,...)$root
	return(soln)
} # nocov end

# call gamma, then take to a power
rgengamma <- function(n,a,d,p) {
	# k, theta are shape, scale
	k <- d / p
	theta <- a ^ p
	retv <- rgamma(n, shape=k, scale=theta)
	retv <- retv ^ (1/p)
	retv
}

# rchisq is currently broken for df=0,ncp=0
# this is a not-recycled chisquare generator
# which irons over the df=0,ncp=0 case.
unbroken_rchisq <- function(n,df,ncp=0) {
	retv <- rchisq(n,df=df,ncp=ncp)
	retv[df==0 & ncp==0] <- 0
	retv
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
