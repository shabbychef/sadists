% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnt.r
\name{dnt}
\alias{dnt}
\alias{ddnt}
\alias{pdnt}
\alias{qdnt}
\alias{rdnt}
\title{The doubly non-central t distribution.}
\usage{
ddnt(x, df, ncp1, ncp2, log = FALSE, order.max=6)

pdnt(q, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)

qdnt(p, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)

rdnt(n, df, ncp1, ncp2)
}
\arguments{
\item{x, q}{vector of quantiles.}

\item{df}{the degrees of freedom for the denominator, \eqn{\nu}{v}.
We do \emph{not} recycle these versus the \code{x,q,p,n}.}

\item{ncp1, ncp2}{the non-centrality parameters for the numerator and denominator,
respectively, \eqn{\mu}{mu} and \eqn{\theta}{theta}
We do \emph{not} recycle these versus the \code{x,q,p,n}.}

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
\code{ddnt} gives the density, \code{pdnt} gives the 
distribution function, \code{qdnt} gives the quantile function, 
and \code{rdnt} generates random deviates.

Invalid arguments will result in return value \code{NaN} with a warning.
}
\description{
Density, distribution function, quantile function and random
generation for the doubly non-central t distribution.
}
\details{
Let \eqn{Z \sim \mathcal{N}\left(\mu,1\right)}{Z ~ N(u,1)} independently
of \eqn{X \sim \chi^2\left(\theta,\nu\right)}{X ~ x^2(theta,v)}. The 
random variable
\deqn{T = \frac{Z}{\sqrt{X/\nu}}}{T = Z / sqrt(X/v)}
takes a \emph{doubly non-central t distribution} with parameters
\eqn{\nu, \mu, \theta}{v, mu, theta}.
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
rvs <- rdnt(128, 20, 1, 1)
dvs <- ddnt(rvs, 20, 1, 1)
pvs.H0 <- pdnt(rvs, 20, 0, 1)
pvs.HA <- pdnt(rvs, 20, 1, 1)
\donttest{
plot(ecdf(pvs.H0))
plot(ecdf(pvs.HA))
}
# compare to singly non-central
dv1 <- ddnt(1, df=10, ncp1=5, ncp2=0, log=FALSE)
dv2 <- dt(1, df=10, ncp=5, log=FALSE)
pv1 <- pdnt(1, df=10, ncp1=5, ncp2=0, log.p=FALSE)
pv11 <- pdnt(1, df=10, ncp1=5, ncp2=0.001, log.p=FALSE)
v2 <- pt(1, df=10, ncp=5, log.p=FALSE)

q1 <- qdnt(pv1, df=10, ncp1=5, ncp2=0, log.p=FALSE)
}
\references{
Krishnan, Marakatha. "Series Representations of the Doubly Noncentral
t-Distribution." Journal of the American Statistical Association 63, 
no. 323 (1968): 1004-1012.
}
\seealso{
t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\keyword{distribution}
