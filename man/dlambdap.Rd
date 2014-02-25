\name{lambdap}
\alias{dlambdap}
\alias{lambdap}
\alias{plambdap}
\alias{qlambdap}
\alias{rlambdap}
\title{The lambda prime distribution.}
\usage{
dlambdap(x, df, t, log = FALSE)

plambdap(q, df, t, lower.tail = TRUE, log.p = FALSE)

qlambdap(p, df, t, lower.tail = TRUE, log.p = FALSE)

rlambdap(n, df, t)
}
\arguments{
  \item{x,q}{vector of quantiles.}

  \item{p}{vector of probabilities.}

  \item{n}{number of observations.}

  \item{df}{the degrees of freedom in the chi square.}

  \item{t}{the scaling parameter on the chi.}

  \item{log}{logical; if TRUE, densities \eqn{f} are given
  as \eqn{\mbox{log}(f)}{log(f)}.}

  \item{log.p}{logical; if TRUE, probabilities p are given
  as \eqn{\mbox{log}(p)}{log(p)}.}

  \item{lower.tail}{logical; if TRUE (default),
  probabilities are \eqn{P[X \le x]}{P[X <= x]}, otherwise,
  \eqn{P[X > x]}{P[X > x]}.}

  \item{x}{vector of blah}

  \item{n}{number of blah}

  \item{...}{arguments passed on to ...}
}
\value{
\code{dlambdap} gives the density, \code{plambdap} gives
the distribution function, \code{qlambdap} gives the
quantile function, and \code{rlambdap} generates random
deviates.

Invalid arguments will result in return value \code{NaN}
with a warning.

\code{dsr} gives the density, \code{psr} gives the
distribution function, \code{qsr} gives the quantile
function, and \code{rsr} generates random deviates.

Invalid arguments will result in return value \code{NaN}
with a warning.
}
\description{
Density, distribution function, quantile function and
random generation for the lambda prime distribution.
}
\details{
Suppose \eqn{y \sim \chi^2\left(\nu\right)}{y ~ x^2(v)},
and \eqn{Z}{Z} is a standard normal. \deqn{T = Z + t
\sqrt{y/\nu}}{T = Z + t sqrt(y/v)} takes a lambda prime
distribution with parameters \eqn{\nu, t}{v, t}. A lambda
prime random variable can be viewed as a confidence
variable on a non-central t because \deqn{t = \frac{Z' +
T}{\sqrt{y/\nu}}}{t = (Z' + T)/sqrt(y/v)}
}
\note{
This is a thin wrapper on the t distribution.
}
\examples{
rvs <- rsr(128, 253*6, 0, 253)
dvs <- dsr(rvs, 253*6, 0, 253)
pvs.H0 <- psr(rvs, 253*6, 0, 253)
pvs.HA <- psr(rvs, 253*6, 1, 253)
\dontrun{
plot(ecdf(pvs.H0))
plot(ecdf(pvs.HA))
}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}

Steven E. Pav \email{steven@cerebellumcapital.com}
}
\references{
Lecoutre, Bruno. "Another look at confidence intervals for
the noncentral t distribution." Journal of Modern Applied
Statistical Methods 6, no. 1 (2007): 107--116.
\url{http://www.univ-rouen.fr/LMRS/Persopage/Lecoutre/telechargements/Lecoutre_Another_look-JMSAM2007_6(1).pdf}

Lecoutre, Bruno. "Two useful distributions for Bayesian
predictive procedures under normal models." Journal of
Statistical Planning and Inference 79 (1999): 93--105.
}
\seealso{
t distribution functions, \code{\link{dt}, \link{pt},
\link{qt}, \link{rt}}, K prime distribution functions,
\code{\link{dkprime}, \link{pkprime}, \link{qkprime},
\link{rkprime}},

t-distribution functions, \code{\link{dt}, \link{pt},
\link{qt}, \link{rt}}
}
\keyword{distribution}
\keyword{io}
\keyword{plotting}

