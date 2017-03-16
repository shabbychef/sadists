#
# dockerfile to CRAN-check with r-dev
#
# docker build --rm -t shabbychef/sadists-crancheck .
#
# docker run -it --rm --volume $(pwd):/srv:rw sadists-crancheck
#
# Created: 2016.01.10
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav
# Comments: Steven E. Pav

#####################################################
# preamble# FOLDUP
FROM shabbychef/crancheck
MAINTAINER Steven E. Pav, shabbychef@gmail.com
# UNFOLD

# rinstall somethings...
RUN /usr/local/bin/install2.r PDQutils hypergeo orthopolynom shiny testthat ggplot2 xtable knitr

#for vim modeline: (do not edit)
# vim:nu:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=Dockerfile:ft=Dockerfile:fo=croql