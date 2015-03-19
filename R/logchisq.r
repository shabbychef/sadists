# /usr/bin/r
#
# Created: 2015.03.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <steven@cerebellumcapital.com>
# Comments: Steven E. Pav

# cumulants of the log of the central chi-square
lc_cumulants <- function(df,order.max=3,orders=c(1:order.max)) {
	kappa <- psigamma(df/2,deriv=orders)
	kappa[1] <- kappa[1] + log(2)
	return(kappa)
}
# moments of the log of the central chi-square
lc_moments <- function(df,order.max=3,orders=c(1:order.max)) {
	kappa <- lc_cumulants(df,orders=orders)
	mu <- PDQutils::cumulant2moment(kappa)
	return(mu)
}
# that not right ...
# moments of the log of the non-central chi-square
lnc_moments <- function(df,ncp=0,order.max=3,orders=c(1:order.max)) {
	allmu <- sapply(0:50,function(iv) {
									(exp(-ncp/2) * ((ncp/2)**(iv)) / factorial(iv)) * lc_moments(df+2*iv,orders=orders)
		},simplify=FALSE)
	mu <- Reduce('+',allmu)
	return(mu)
}
# cumulants of the log of the non-central chi-square
lnc_cumulants <- function(df,ncp=0,order.max=3,orders=c(1:order.max)) {
	mu <- lnc_moments(df,ncp,orders=orders)
	kappa <- PDQutils::moment2cumulant(mu)
	return(kappa)
}






#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
