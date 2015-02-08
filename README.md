

# sadists

Some Auxiliary Distributions apparently not available in R.

-- Steven E. Pav, shabbychef@gmail.com

## Installation

This package may be installed from CRAN; the latest version may be
found on [github](https://www.github.com/shabbychef/sadists "sadists")
via devtools:


```r
if (require(devtools)) {
    # latest greatest
    install_github(repo = "sadists", username = "shabbychef", 
        ref = "master")
}
```

## Doubly non-central t distribution

The [doubly non-central t distribution](https://en.wikipedia.org/wiki/Doubly_noncentral_t-distribution)
generalizes the t distribution to the case where the denominatory
chi-square is non-central.


```r
k <- 5
mu <- 1
theta <- 2
rvs <- rdnt(1000, k, mu, theta)
```

```
## Error in eval(expr, envir, enclos): could not find function "rdnt"
```

```r
pvs <- pdnt(rvs, k, mu, theta)
```

```
## Error in eval(expr, envir, enclos): could not find function "pdnt"
```

```r
qvs <- qdnt(pvs, k, mu, theta)
```

```
## Error in eval(expr, envir, enclos): could not find function "qdnt"
```

```r
plot(ecdf(pvs))
```

```
## Error in sort(x): object 'pvs' not found
```




# Basic Usage

