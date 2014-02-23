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

# quantile function helper; 
# call uniroot in a neighborhood of x0
# the function fnc is non-decreasing
uniroot_helper <- function(fnc, x0=0, f0=NULL, tol=.Machine$double.eps^0.25,
													 maxiter=1000, ...) {
	if (is.null(f0)) 
		f0 = fnc(x0,...)
	if (f0 == 0)
		return(x0)
	if (f0 > 0) {
		ub <- x0
		fub <- f0
		lb <- x0 - 0.1 * max(0.01,abs(x0))
		flb <- fnc(lb,...)
		# check for NA, infs, etc.
		while (flb > 0) {
			# drag ub down too
			ub <- lb
			fub <- flb
			lb <- q0 - 2 * (q0 - lb)
			flb <- fnc(lb,...)
		}
	} else {
		lb <- q0
		flb <- v0
		ub <- q0 + 0.1 * max(0.01,abs(q0))
		fub <- fnc(ub,...)
# 2FIX: beware infs!
		while (fub < 0) {
			# drag lb up too
			lb <- ub
			flb <- fub
			ub <- q0 + 2 * (ub - q0)
			fub <- fnc(ub,...)
		}
	}
	soln <- uniroot(zerf,c(lb,ub),f.lower=flb,f.upper=fub,tol=tol,maxiter=maxiter,...)$root
	return(soln)
}


#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
