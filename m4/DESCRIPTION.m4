dnl divert here just means the output from basedefs does not appear.
divert(-1)
include(basedefs.m4)
divert(0)dnl
Package: PKG_NAME()
Maintainer: Steven E. Pav <shabbychef@gmail.com>
Authors@R: c(person(c("Steven", "E."), "Pav", 
    role=c("aut","cre"),
    email="shabbychef@gmail.com",
    comment = c(ORCID = "0000-0002-4197-6195")))
Version: VERSION()
Date: DATE()
License: LGPL-3
Title: Some Additional Distributions
BugReports: https://github.com/shabbychef/PKG_NAME()/issues
Description: Provides the density, distribution, quantile and generation functions of some obscure probability 
    distributions, including the doubly non-central t, F, Beta, and Eta distributions; 
    the lambda-prime and K-prime; the upsilon distribution; the (weighted) sum of 
    non-central chi-squares to a power; the (weighted) sum of log non-central chi-squares;
    the product of non-central chi-squares to powers; the product of doubly non-central
    F variables; the product of independent normals.
Depends: 
    R (>= 3.0.2)
Imports:
    PDQutils (>= 0.1.1),
    hypergeo,
    orthopolynom
Suggests: 
    SharpeR,
    shiny,
    testthat, 
    ggplot2, 
    xtable,
    formatR,
    knitr
URL: https://github.com/shabbychef/PKG_NAME()
VignetteBuilder: knitr
Collate:
m4_R_FILES()
dnl vim:ts=4:sw=2:tw=79:syn=m4:ft=m4:et
