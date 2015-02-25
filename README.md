

# sadists

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
pvs <- pdnt(rvs, k, mu, theta)
qvs <- qdnt(pvs, k, mu, theta)
```

```
## Error: object 'q0' not found
```

```r
plot(ecdf(pvs))
```

![plot of chunk dnt](github_extra/figure/dnt-1.png) 

## K-prime distribution

## Lambda prime distribution

## Upsilon distribution



# Basic Usage

