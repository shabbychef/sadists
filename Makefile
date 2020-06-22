######################
# 
# Created: 2016.03.26
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
######################

############### FLAGS ###############

PKG_NAME 					:= sadists
VMAJOR 						 = 0
VMINOR 						 = 2
VPATCH  					 = 3
VDEV 							 = .3000

undefine RPKG_USES_RCPP

include ./rpkg_make/Makefile

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=.tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
