

# sadists

[![Build Status](https://travis-ci.org/shabbychef/sadists.png)](https://travis-ci.org/shabbychef/sadists)

Some Additional Distributions apparently not available in R.

-- Steven E. Pav, shabbychef@gmail.com

## Installation

This package may be installed from CRAN; the latest version may be
found on [github](https://www.github.com/shabbychef/sadists "sadists")
via devtools, or installed via [drat](https://github.com/eddelbuettel/drat "drat"):


```r
if (require(devtools)) {
    # latest greatest
    install_github("shabbychef/sadists")
}
# via drat:
if (require(drat)) {
    drat:::add("shabbychef")
    install.packages("sadists")
}
```

## Doubly non-central t distribution

The [doubly non-central t distribution](https://en.wikipedia.org/wiki/Doubly_noncentral_t-distribution)
generalizes the t distribution to the case where the denominator chi-square is non-central.


```r
require(sadists)
k <- 5
mu <- 1
theta <- 2
rvs <- rdnt(1000, k, mu, theta)
# pvs <- pdnt(rvs, k, mu, theta) qvs <- qdnt(pvs,
# k, mu, theta) plot(ecdf(pvs))
```

## Doubly non-central F distribution

The doubly non-central F distribution generalizes the F distribution to the case where the denominator
chi-square is non-central. It has not yet been implemented.

## Sum of (non-central) Chi-squares

The weighted sum of chi-squares is not a commonly seen random variable. However, its cumulants can be
easily computed, so its 'PDQ' functions can easily be computed. Moreover, its distribution and quantile
functions can be used in computation of those of the doubly non-central F.


```r
require(sadists)
wts <- c(1, -3, 4)
df <- c(100, 20, 10)
ncp <- c(5, 3, 1)
rvs <- rsumchisq(128, wts, df, ncp)
dvs <- dsumchisq(rvs, wts, df, ncp)
qvs <- psumchisq(rvs, wts, df, ncp)
pvs <- qsumchisq(ppoints(length(rvs)), wts, df, ncp)
plot(ecdf(qvs))
```

![plot of chunk schisq](github_extra/figure/schisq-1.png) 

## K-prime distribution

## Lambda prime distribution

## Upsilon distribution



# Basic Usage

