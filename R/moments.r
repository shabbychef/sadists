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
# Author: Steven E. Pav
# Comments: Steven E. Pav

# utilities#FOLDUP
# central moments to standardized moments
central2std <- function(mu.cent) {
	#std <- sqrt(mu.cent[3])
	mu.std <- mu.cent / (mu.cent[3] ^ ((0:(length(mu.cent)-1))/2))
	return(mu.std)
}
#UNFOLD

# some moments of standard distributions:#FOLDUP

# compute the 1 through order.max raw, uncentered moment
# of the normal distribution with given mean and standard
# deviation
norm.moms <- function(mu=0,sigma=1,order.max=3) {
	retval <- rep(1,order.max)
	hermi <- orthopolynom::hermite.he.polynomials(order.max, normalized=FALSE)
	for (iii in c(1:order.max)) {
		cvals <- abs(coefficients(hermi[[iii+1]]))
		lvals <- mu^(0:iii) * sigma^(iii - (0:iii))
		retval[iii] <- sum(cvals * lvals)
	}
	return(retval)
}

# compute the 1 through order.max raw, uncentered moment
# of the (central) chi distribution with df d.f.
chi.moms <- function(df,order.max=3,orders=1:order.max,log=FALSE) {
	retval <- chisq.moms(df=df,orders=orders/2,log=log)
	return(retval)
}

# compute the 1 through order.max raw, uncentered moment
# of the (central) chi-square distribution with df d.f.
chisq.moms <- function(df,order.max=3,orders=1:order.max,log=FALSE) {
	if (log) {
		retval <- orders * log(2) + lgamrat(orders + (df/2),df/2)
	} else {
		retval <- (2^(orders)) * gamrat(orders + (df/2),df/2)
	}
	return(retval)
}

# something like a nakagami, but really a scaled chi
schi.moms <- function(df,scal=1,order.max=3) {
	stopifnot(df > 0)
	orders <- 1:order.max
	if (is.infinite(df)) {
		retval <- scal ^ orders
	} else {
		retval <- chi.moms(df=df,orders=orders,log=TRUE)
		retval <- retval + orders * (log(abs(scale)) - 0.5 * log(df))
		retval <- exp(retval) * sign(scale)^orders
	}
	return(retval)
}


#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
