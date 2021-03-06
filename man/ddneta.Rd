% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dneta.r
\name{dneta}
\alias{dneta}
\alias{ddneta}
\alias{pdneta}
\alias{qdneta}
\alias{rdneta}
\title{The doubly non-central Eta distribution.}
\usage{
ddneta(x, df, ncp1, ncp2, log = FALSE, order.max=6)

pdneta(q, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)

qdneta(p, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)

rdneta(n, df, ncp1, ncp2)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{df}{the degrees of freedom for the denominator chi square.
We do \emph{not} recycle this versus the \code{x,q,p,n}.}

\item{ncp1, ncp2}{the non-centrality parameters for the numerator and denominator.
We do \emph{not} recycle these versus the \code{x,q,p,n}.
Note that the sign of \code{ncp1} is important, while
\code{ncp2} must be non-negative.}

\item{log}{logical; if TRUE, densities \eqn{f} are given 
as \eqn{\mbox{log}(f)}{log(f)}.}

\item{order.max}{the order to use in the approximate density,
distribution, and quantile computations, via the Gram-Charlier,
Edeworth, or Cornish-Fisher expansion.}

\item{p}{vector of probabilities.}

\item{n}{number of observations.}

\item{log.p}{logical; if TRUE, probabilities p are given 
as \eqn{\mbox{log}(p)}{log(p)}.}

\item{lower.tail}{logical; if TRUE (default), probabilities are
\eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
\code{ddneta} gives the density, \code{pdneta} gives the 
distribution function, \code{qdneta} gives the quantile function, 
and \code{rdneta} generates random deviates.

Invalid arguments will result in return value \code{NaN} with a warning.
}
\description{
Density, distribution function, quantile function and random
generation for the doubly non-central Eta distribution.
}
\details{
Suppose \eqn{Z}{Z} is a normal with mean \eqn{\delta_1}{delta_1},
and standard deviation 1, independent of 
\eqn{X \sim \chi^2\left(\delta_2,\nu_2\right)}{X ~ X^2(delta_2,v_2)},
a non-central chi-square with \eqn{\nu_2}{v_2} degrees of freedom
and non-centrality parameter \eqn{\delta_2}{delta_2}. Then
\deqn{Y = \frac{Z}{\sqrt{Z^2 + X}}}{Y = Z/sqrt(Z^2 + X)}
takes a doubly non-central Eta distribution with 
\eqn{\nu_2}{v_2} degrees of freedom and non-centrality parameters
\eqn{\delta_1,\delta_2}{delta_1,delta_2}. The \emph{square} of
a doubly non-central Eta is a doubly non-central Beta variate.
}
\note{
The PDF, CDF, and quantile function are approximated, via
the Edgeworth or Cornish Fisher approximations, which may
not be terribly accurate in the tails of the distribution.
You are warned.

The distribution parameters are \emph{not} recycled
with respect to the \code{x, p, q} or \code{n} parameters,
for, respectively, the density, distribution, quantile
and generation functions. This is for simplicity of
implementation and performance. It is, however, in contrast
to the usual R idiom for dpqr functions.
}
\examples{
rv <- rdneta(500, df=100,ncp1=1.5,ncp2=12)
d1 <- ddneta(rv, df=100,ncp1=1.5,ncp2=12)
\donttest{
plot(rv,d1)
}
p1 <- ddneta(rv, df=100,ncp1=1.5,ncp2=12)
# should be nearly uniform:
\donttest{
plot(ecdf(p1))
}
q1 <- qdneta(ppoints(length(rv)), df=100,ncp1=1.5,ncp2=12)
\donttest{
qqplot(x=rv,y=q1)
}
}
\seealso{
(doubly non-central) t distribution functions, 
\code{\link{ddnt}, \link{pdnt}, \link{qdnt}, \link{rdnt}}.

(doubly non-central) Beta distribution functions, 
\code{\link{ddnbeta}, \link{pdnbeta}, \link{qdnbeta}, \link{rdnbeta}}.
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\keyword{distribution}
