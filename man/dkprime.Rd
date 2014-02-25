\name{kprime}
\alias{dkprime}
\alias{kprime}
\alias{pkprime}
\alias{qkprime}
\alias{rkprime}
\title{The K prime distribution.}
\usage{
dkprime(x, v1, v2, a, b = 1, log = FALSE)

pkprime(q, v1, v2, a, b = 1, lower.tail = TRUE, log.p = FALSE)

qkprime(p, v1, v2, a, b = 1, lower.tail = TRUE, log.p = FALSE)

rkprime(n, v1, v2, a, b = 1)
}
\arguments{
  \item{x,q}{vector of quantiles.}

  \item{p}{vector of probabilities.}

  \item{n}{number of observations.}

  \item{v1}{the degrees of freedom in the numerator
  chisquare. When (positive) infinite, we recover a
  non-central t distribution with \code{v2} degrees of
  freedom and non-centrality parameter \code{a}, scaled by
  \code{b}.}

  \item{v2}{the degrees of freedom in the denominator
  chisquare. When equal to infinity, we recover the Lambda
  prime distribution.}

  \item{a}{the non-centrality scaling parameter. When equal
  to zero, we recover the (central) t distribution.}

  \item{b}{the scaling parameter.}

  \item{log}{logical; if TRUE, densities \eqn{f} are given
  as \eqn{\mbox{log}(f)}{log(f)}.}

  \item{log.p}{logical; if TRUE, probabilities p are given
  as \eqn{\mbox{log}(p)}{log(p)}.}

  \item{lower.tail}{logical; if TRUE (default),
  probabilities are \eqn{P[X \le x]}{P[X <= x]}, otherwise,
  \eqn{P[X > x]}{P[X > x]}.}
}
\value{
\code{dkprime} gives the density, \code{pkprime} gives the
distribution function, \code{qkprime} gives the quantile
function, and \code{rkprime} generates random deviates.

Invalid arguments will result in return value \code{NaN}
with a warning.
}
\description{
Density, distribution function, quantile function and
random generation for the K prime distribution.
}
\details{
Suppose \eqn{y \sim \chi^2\left(\nu_1\right)}{y ~ x^2(v1)},
and \eqn{x \sim t \left(\nu_2, a\sqrt{y/\nu_1}/b\right)}{x
~ t(v2,(a/b) sqrt(y/v1))}. Then the random variable \deqn{T
= b x}{T = b x} takes a K prime distribution with
parameters \eqn{\nu_1, \nu_2, a, b}{v1, v2, a, b}. In
Lecoutre's terminology, \eqn{T \sim K'_{\nu_1,
\nu_2}\left(a, b\right)}{T ~ K'_v1,v2(a,b)}

Equivalently, we can think of \deqn{T = \frac{b Z + a
\sqrt{\chi^2_{\nu_1} / \nu_1}}{\sqrt{\chi^2_{\nu_2} /
\nu_2}}}{T = (bZ + a sqrt(chi2_v1/v1)) / sqrt(chi2_v2/v2)}
where \eqn{Z} is a standard normal, and the normal and the
(central) chi-squares are independent of each other. When
\eqn{a=0}{a=0} we recover a central t distribution; when
\eqn{\nu_1=\infty}{v1=inf} we recover a rescaled
non-central t distribution; when \eqn{b=0}{b=0}, we get a
rescaled square root of a central F distribution; when
\eqn{\nu_2=\infty}{v2=inf}, we recover a Lambda prime
distribution.
}
\examples{
d1 <- dkprime(1, 50, 20, a=0.01)
d2 <- dkprime(1, 50, 20, a=0.0001)
d3 <- dkprime(1, 50, 20, a=0)
d4 <- dkprime(1, 10000, 20, a=1)
d5 <- dkprime(1, Inf, 20, a=1)

\dontrun{
avals <- 10^(seq(1,7,length.out=101))
plot(log10(avals),dvals)
}
}
\author{
Steven E. Pav \email{shabbychef@gmail.com}
}
\references{
Lecoutre, Bruno. "Two Useful distributions for Bayesian
predictive procedures under normal models." Journal of
Statistical Planning and Inference 79, no. 1 (1999):
93-105.

Poitevineau, Jacques, and Lecoutre, Bruno. "Implementing
Bayesian predictive procedures: The K-prime and K-square
distributions." Computational Statistics and Data Analysis
54, no. 3 (2010): 724-731.
\url{http://arxiv.org/abs/1003.4890v1}
}
\seealso{
t distribution functions, \code{\link{dt}, \link{pt},
\link{qt}, \link{rt}}, K square distribution functions,
\code{\link{dksquare}, \link{pksquare}, \link{qksquare},
\link{rksquare}}, lambda prime distribution functions,
\code{\link{dlambdap}, \link{plambdap}, \link{qlambdap},
\link{rlambdap}},
}
\keyword{distribution}

