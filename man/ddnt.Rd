\name{dnt}
\alias{ddnt}
\alias{dnt}
\alias{pdnt}
\alias{qdnt}
\alias{rdnt}
\title{The doubly non-central t distribution.}
\usage{
ddnt(x, k, mu, theta, log = FALSE)

pdnt(q, k, mu, theta, lower.tail = TRUE, log.p = FALSE)

qdnt(p, k, mu, theta, lower.tail = TRUE, log.p = FALSE)

rdnt(n, k, mu, theta)
}
\arguments{
  \item{x,q}{vector of quantiles.}

  \item{p}{vector of probabilities.}

  \item{n}{number of observations.}

  \item{k}{the degrees of freedom in the denominator.}

  \item{mu}{the numerator non-centrality parameter,
  \eqn{\mu}{mu}.}

  \item{theta}{the denominator non-centrality parameter,
  \eqn{\theta}{theta}. When equal to zero, we recover the
  singly non-central t distribution.}

  \item{log}{logical; if TRUE, densities \eqn{f} are given
  as \eqn{\mbox{log}(f)}{log(f)}.}

  \item{log.p}{logical; if TRUE, probabilities p are given
  as \eqn{\mbox{log}(p)}{log(p)}.}

  \item{lower.tail}{logical; if TRUE (default),
  probabilities are \eqn{P[X \le x]}{P[X <= x]}, otherwise,
  \eqn{P[X > x]}{P[X > x]}.}
}
\value{
\code{ddnt} gives the density, \code{pdnt} gives the
distribution function, \code{qdnt} gives the quantile
function, and \code{rdnt} generates random deviates.

Invalid arguments will result in return value \code{NaN}
with a warning.
}
\description{
Density, distribution function, quantile function and
random generation for the doubly non-central t
distribution.
}
\details{
Let \eqn{X \sim \mathcal{N}\left(\mu,1\right)}{X ~ N(u,1)}
independently of \eqn{Y \sim \chi^2\left(k,\theta\right)}{Y
~ x^2(k,theta)}. The random variable \deqn{T =
\frac{X}{\sqrt{Y/k}}}{T = X / sqrt(Y/k)} takes a
\emph{doubly non-central t distribution} with parameters
\eqn{k, \mu, \theta}{k, mu, theta}.
}
\examples{
rvs <- rdnt(128, 20, 1, 1)
dvs <- ddnt(rvs, 20, 1, 1)
pvs.H0 <- pdnt(rvs, 20, 0, 1)
pvs.HA <- pdnt(rvs, 20, 1, 1)
\dontrun{
plot(ecdf(pvs.H0))
plot(ecdf(pvs.HA))
}
# compare to singly non-central
dv1 <- ddnt(1, k=10, mu=5, theta=0, log=FALSE)
dv2 <- dt(1, df=10, ncp=5, log=FALSE)
pv1 <- pdnt(1, k=10, mu=5, theta=0, log.p=FALSE)
pv11 <- pdnt(1, k=10, mu=5, theta=0.001, log.p=FALSE)
v2 <- pt(1, df=10, ncp=5, log.p=FALSE)

q1 <- qdnt(pv1, k=10, mu=5, theta=0, log.p=FALSE)
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\references{
Krishnan, Marakatha. "Series Representations of the Doubly
Noncentral t-Distribution." Journal of the American
Statistical Association 63, no. 323 (1968): 1004-1012.
}
\seealso{
t distribution functions, \code{\link{dt}, \link{pt},
\link{qt}, \link{rt}}
}
\keyword{distribution}

