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

# Created: 2015.02.23
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

# compute the 1 through order.max raw, cumulants 
# of the normal distribution with given mean and sd.
norm.cumuls <- function(mean=0,sd=1,order.max=3) {
	retval <- rep(0,order.max)
	if (order.max > 0)
		retval[1] <- mean
	if (order.max > 1)
		retval[2] <- sd
	return(retval)
}

# compute the 1 through order.max raw cumulants
# of the (central) chi distribution with df d.f.
chi.cumults <- function(df,order.max=3) {
	retval <- moment2cumulant(chi.moms(df,order.max=order.max))
	return(retval)
}

# compute the 1 through order.max raw, cumulants 
# of the (non-central) chi-square distribution with df d.f.
# and noncentrality parameter ncp
chisq.cumuls <- function(df,ncp=0,order.max=3) {
	jvals <- 0:(order.max-1)
	retval <- (2^(jvals)) * factorial(jvals) * (df + ncp * (jvals+1))
	return(retval)
}

# something like a nakagami, but really a scaled chi
schi.cumuls <- function(df,scal=1,order.max=3) {
	retval <- moment2cumulant(schi.moms(df,scal=scal,order.max=order.max))
	return(retval)
}



#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
