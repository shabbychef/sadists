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
FROM rocker/r-base:latest
MAINTAINER Steven E. Pav, shabbychef@gmail.com
# UNFOLD

ENV DOCKERFILE_REFRESHED_AT 2020.06.23
# see http://crosbymichael.com/dockerfile-best-practices.html
#RUN echo "deb http://archive.ubuntu.com/ubuntu precise main universe" > /etc/apt/sources.list

RUN (apt-get clean -y ; \
 apt-get update -y -qq; \
 apt-get update --fix-missing );

# from https://superuser.com/a/1557981/60965
RUN (DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true apt-get install -q -y --no-install-recommends \ 
  gcc-9-base libgcc-9-dev libc6-dev ; )

#RUN (apt-get dist-upgrade -y ; \
#apt-get update -qq ; \
RUN (DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true apt-get install -q -y --no-install-recommends \ 
  libgs9 texlive-base texlive-binaries libcupsimage2 libcups2 file curl wget qpdf pandoc ghostscript \
  texlive-latex-extra texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texinfo \
  libmagick++-dev libmariadbclient-dev libpq-dev liblapack-dev libblas-dev libssl-dev libcurl4-openssl-dev libgit2-dev libxml2-dev libcairo2-dev; \
  sync ; \
  mkdir -p /usr/local/lib/R/site-library ; \
  chmod -R 777 /usr/local/lib/R/site-library ; \
  sync )

RUN (DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true apt-get install -q -y --no-install-recommends \ 
  cargo mesa-common-dev ; \
  sync )

# rinstall somethings...
RUN (install2.r --error --skipinstalled --deps TRUE \
  devtools roxygen2 testthat curl \
  knitr formatR codetools; )

# package dependencies
RUN (install2.r --error --skipinstalled --deps TRUE \
	PDQutils hypergeo orthopolynom SharpeR \
	shiny testthat ggplot2 xtable ; )

# already in there, I think. weird.
# RUN groupadd -g 1000 spav && useradd -g spav -u 1000 spav;
# USER spav
# docker has id 1000 in rocker/r-base; so use that. ugh.
USER docker

WORKDIR /srv

#####################################################
# entry and cmd# FOLDUP
# these are the default, but remind you that you might want to use /usr/bin/R instead?
# always use array syntax:
ENTRYPOINT ["R","CMD","check","--as-cran","--output=/tmp"]

# ENTRYPOINT and CMD are better together:
CMD ["/srv/*.tar.gz"]
# UNFOLD

#for vim modeline: (do not edit)
# vim:nu:et:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=Dockerfile:ft=Dockerfile:fo=croql
