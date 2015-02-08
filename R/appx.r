# /usr/bin/r
#
# Created: 2015.02.07
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav


# for the Hermite Polynomials
require(orthopolynom)
require(moments)

# utilities or dealing with moments and cumulants

#' @title Convert moments to raw cumulants.
#'
#' @description 
#'
#' Conversion of a vector of moments to raw cumulants.
#'
#' @details
#'
#' The 'raw' cumulants \eqn{\kappa_i}{kappa_i} are connected
#' to the 'raw' (uncentered) moments, \eqn{\mu_i'}{mu'_i} via
#' the equation
#' \deqn{\kappa_n = \mu_n' - \sum_{m=1}^{n-1} {n-1 \choose m-1} \kappa_m \mu_{n-m}'}
#'
#' Note that this formula also works for central moments, assuming
#' the distribution has been normalized to zero mean.
#'
#' @usage
#'
#' moment2cumulant(moms)
#'
#' @param moms a vector of the moments. The first element is the first moment.
#' If centered moments are given, the first cumulant shall be zero.
#' @return a vector of the cumulants.
#'
#' @keywords distribution 
#' @seealso \code{\link{cumulant2moment}}
#' @export 
#'
#' @examples 
#' # normal distribution, mean 0, variance 1
#' n.cum <- moment2cumulant(c(0,1,0,3,0,15))
#' # normal distribution, mean 1, variance 1
#' n.cum <- moment2cumulant(c(1,2,4,10,26))
#' # exponential distribution
#' lambda <- 0.7
#' n <- 1:6
#' e.cum <- moment2cumulant(factorial(n) / (lambda^n))
#' @template etc
moment2cumulant <- function(moms) {
	kappa <- moms
	if (length(kappa) > 1) {
		for (nnn in 2:length(kappa)) {
			mmm <- 1:(nnn-1)
			kappa[nnn] <- moms[nnn] - sum(choose(nnn-1,mmm-1) * kappa[mmm] * moms[nnn-mmm])
		}
	}
	return(kappa)
}

#' @title Convert raw cumulants to moments.
#'
#' @description 
#'
#' Conversion of a vector of raw cumulatnts to moments.
#'
#' @details
#'
#' The 'raw' cumulants \eqn{\kappa_i}{kappa_i} are connected
#' to the 'raw' (uncentered) moments, \eqn{\mu_i'}{mu'_i} via
#' the equation
#' \deqn{\mu_n' = \kappa_n + \sum_{m=1}^{n-1} {n-1 \choose m-1} \kappa_m \mu_{n-m}'}
#'
#' @usage
#'
#' cumulant2moment(kappa)
#'
#' @param kappa a vector of the raw cumulants. The first element is the first cumulant,
#' which is also the first moment.
#' @return a vector of the raw moments.
#'
#' @keywords distribution 
#' @seealso \code{\link{moment2cumulant}}
#' @export 
#'
#' @examples 
#' # normal distribution, mean 0, variance 1
#' n.mom <- cumulant2moment(c(0,1,0,0,0,0))
#' # normal distribution, mean 1, variance 1
#' n.mom <- cumulant2moment(c(1,1,0,0,0,0))
#' @template etc
cumulant2moment <- function(kappa) {
	moms <- kappa
	if (length(moms) > 1) {
		for (nnn in 2:length(kappa)) {
			mmm <- 1:(nnn-1)
			moms[nnn] <- kappa[nnn] + sum(choose(nnn-1,mmm-1) * kappa[mmm] * moms[nnn-mmm])
		}
	}
	return(moms)
}

# central moments to standardized moments
central2std <- function(mu.cent) {
	#std <- sqrt(mu.cent[3])
	mu.std <- mu.cent / (mu.cent[3] ^ ((0:(length(mu.cent)-1))/2))
	return(mu.std)
}

#' @title Higher order Cornish Fisher approximation.
#'
#' @description 
#'
#' Cornish Fisher approximate quantiles on a standardized random variable.
#'
#' @details
#'
#' This implements Algorithm AS269 of Lee and Lin for computing higher
#' order Cornish Fisher approximates of the quantile of a standardized
#' random variable. The Cornish Fisher approximation is the Legendre
#' inversion of the Edgeworth expansion of a distribution, but ordered
#' in a way that is convenient when used on the mean of a number of
#' independent draws of a random variable. 
#'
#' Suppose \eqn{x_1, x_2, \ldots, x_n}{x_1, x_2, ..., x_n} are \eqn{n} independent 
#' draws from some probability distribution. 
#' Letting 
#' \deqn{X = \frac{1}{\sqrt{n}} \sum_{1 \le i \le n} x_i,}{X = (x_1 + x_2 + ... x_n) / sqrt(n),}
#' the Central Limit Theorem assures us that, assuming finite variance, 
#' \deqn{X \rightsquigarrow \mathcal{N}(\sqrt{n}\mu, \sigma),}{X ~~ N(sqrt(n) mu, sigma),}
#' with convergence in \eqn{n}.
#'
#' The Cornish Fisher approximation gives a more detailed picture of the
#' quantiles of \eqn{X}{X}, one that is arranged in decreasing powers of
#' \eqn{\sqrt{n}}{sqrt(n)}. The quantile function is the function \eqn{q(p)}{q(p)} 
#' such that \eqn{P\left(X \le q(p)\right) = q(p)}{P(x <= q(p)) = p}. The
#' Cornish Fisher expansion is 
#' \deqn{q(p) = \sqrt{n}\mu + \sigma \left(z + \sum_{3 \le j} c_j f(z)\right),}{q(p) = sqrt{n}mu + sigma (z + sum_{3 <= j} c_j f(z)),}
#' where \eqn{z = \Phi^{-1}(p)}{z = qnorm(p)}, and \eqn{c_j}{c_j} involves
#' standardized cumulants of the distribution of \eqn{x_i}{x_i} of order
#' \eqn{j} and higher. Moreover, the \eqn{c_j}{c_j} feature decreasing powers
#' of \eqn{\sqrt{n}}{sqrt(n)}, giving some justification for truncation.
#' When \eqn{n=1}{n=1}, however, the ordering is somewhat arbitrary.
#'
#' @usage
#'
#' AS269(z,cumul,order.max=NULL,all.ords=FALSE)
#'
#' @param z the quantiles of the normal distribution. an atomic vector.
#' @param cumul the standardized cumulants of order 3, 4, ..., k. an atomic
#' vector.
#' @param order.max the maximum order approximation, must be greater than
#' \code{length(cumul)+2}.
#' We assume the cumulants have been adjusted to reflect that the random
#' variable has unit variance ('standardized cumulants')
#' @param all.ords a logical value. If \code{TRUE}, then results are returned
#' as a matrix, with a column for each order of the approximation. Otherwise
#' the results are a matrix with a single column of the highest order
#' approximation.
#' @return A matrix, which is, depending on \code{all.ords}, either with one column per 
#' order of the approximation, or a single column giving the maximum order
#' approximation. There is one row per value in \code{z}.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @keywords distribution 
#' @seealso \code{\link{qapx.cf}}
#' @export 
#'
#' @template ref-AS269
#' @template ref-Jaschke
#'
#' @examples 
#' foo <- AS269(seq(-2,2,0.01),c(0,2,0,4))
#' # test with the normal distribution:
#' s.cumul <- c(0,0,0,0,0,0,0,0,0)
#' pv <- seq(0.001,0.999,0.001)
#' zv <- qnorm(pv)
#' apq <- AS269(zv,s.cumul,all.ords=FALSE)
#' err <- zv - apq
#'
#' # test with the exponential distribution
#' rate <- 0.7
#' n <- 18
#' # these are 'raw' cumulants'
#' cumul <- (rate ^ -(1:n)) * factorial(0:(n-1))
#' # standardize and chop
#' s.cumul <- cumul[3:length(cumul)] / (cumul[2]^((3:length(cumul))/2))
#' pv <- seq(0.001,0.999,0.001)
#' zv <- qnorm(pv)
#' apq <- cumul[1] + sqrt(cumul[2]) * AS269(zv,s.cumul,all.ords=TRUE)
#' truq <- qexp(pv, rate=rate)
#' err <- truq - apq
#' colSums(abs(err))
#'
#' # an example from Wikipedia page on CF, \url{https://en.wikipedia.org/wiki/Cornish%E2%80%93Fisher_expansion}
#' s.cumul <- c(5,2)
#' apq <- 10 + sqrt(25) * AS269(qnorm(0.95),s.cumul,all.ords=TRUE)
#'
#' @template etc
AS269 <- function(z,cumul,order.max=NULL,all.ords=FALSE) {#FOLDUP
	if (is.null(order.max))
		order.max <- length(cumul) + 2
	if (order.max <= 2)
		return(z)
	stopifnot(order.max <= length(cumul)+2)

	# I use lower case letters for variables which are scalar
	# in z (or, rather, 'broadcast' across all z), and 
	# UPPER CASE for those which are the same size as 
	# z in one dimension or another. 
	# Otherwise, I have tried to make the variables match
	# the published AS269.
	nord <- order.max - 2

	# pre-line 10
	jidx <- 1:nord
	a <- ((-1)^jidx) * cumul / ((jidx+1) * (jidx+2))
	
	# line 10
	# prealloc H 
	n <- length(z)
	H <- matrix(0,nrow=n,ncol=3*nord)
	H[,1] <- - z
	H[,2] <- z * z - 1
	for (jj in 3:(dim(H)[2]))
		H[,jj] <- -(z * H[,jj-1] + (jj-1) * H[,jj-2])

	# line 20
	P <- matrix(0,nrow=n,ncol=3*nord*(nord+1)/2)

	# line 30
	D <- matrix(0,nrow=n,ncol=3*nord)
	DEL <- matrix(0,nrow=n,ncol=ifelse(all.ords,nord,1))

	D[,1] <- -a[1] * H[,2]
	DEL[,1] <- D[,1]
	P[,1] <- D[,1]
	P[,3] <- a[1]
	ja <- 0
	fac <- 1
	
	# main loop
	if (nord > 1) {
		for (jj in 2:nord) {#FOLDUP
			fac <- fac * jj
			ja <- ja + 3 * (jj-1)
			jb <- ja
			bc <- 1
			# calculate coefficients of Hermite polynomitals
			for (kk in 1:(jj-1)) {
				DD <- bc * D[,kk]
				aa <- bc * a[kk]
				jb <- jb - 3 * (jj - kk)
				for (ll in (1:(3*(jj-kk)))) {
					jbl <- jb + ll
					jal <- ja + ll
					P[,jal+1] <- P[,jal+1] + DD * P[,jbl]
					P[,jal+kk+2] <- P[,jal+kk+2] + aa * P[,jbl]
				}  # line 40
				bc <- bc * (jj - kk) / kk
			}  # line 50
			P[,ja + jj + 2] <- P[,ja + jj + 2] + a[jj]
			# calculate the adjustments
			D[,jj] <- 0
			for (ll in (2:(3*jj))) 
				D[,jj] <- D[,jj] - P[,ja + ll] * H[,ll-1]
			# line 60
			P[,ja+1] <- D[,jj]
			if (all.ords) 
				DEL[,jj] <- D[,jj] / fac
			else
				DEL <- DEL + (D[,jj] / fac)
		}  # line 70#UNFOLD
	}

	# the quantile approximations are then 
	# z + columns of DEL
	if (all.ords) 
		retval <- (z + t(apply(t(DEL),2,cumsum)))
	else
		retval <- z + DEL

	return(retval)
}#UNFOLD

#' @title Approximate density and distribution via Gram-Charlier A expansion.
#'
#' @description 
#'
#' Approximate the probability density or cumulative distribution function of a distribution via its raw moments.
#'
#' @details
#'
#' Given the raw moments of a probability distribution, we approximate the probability 
#' density, or the cumulative distribution function, via a Gram-Charlier A expansion on the 
#' standardized distribution.
#'
#' Suppose \eqn{f(x)}{f(x)} is the probability density of some random
#' variable, and let \eqn{F(x)}{F(x)} be the cumulative distribution function.
#' Let \eqn{He_j(x)}{He_j(x)} be the \eqn{j}{j}th probabilist's Hermite
#' polynomial. These polynomials form an orthogonal basis, with respect to the
#' function \eqn{w(x)}{w(x)} of the Hilbert space of functions which are square
#' \eqn{w}{w}-weighted integrable. The weighting function is 
#' \eqn{w(x) = e^{-x^2/2} = \sqrt{2\pi}\phi(x)}{w(x) = e^{-x^2/2} = sqrt(2pi) phi(x)}.
#' The orthogonality relationship is
#' \deqn{\int_{-\infty}^{\infty} He_i(x) He_j(x) w(x) \mathrm{d}x = \sqrt{2\pi} j! \delta_{ij}.}{integral_-inf^inf He_i(x) He_j(x) w(x) dx = sqrt(2pi)j!dirac_ij.}
#'
#' Expanding the density \eqn{f(x)}{f(x)} in terms of these polynomials in the
#' usual way (abusing orthogonality) one has
#' \deqn{f(x) = \sum_{0\le j} \frac{He_j(x)}{j!} \phi(x) \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{f(x) = sum_{0 <= j} (He_j(x)/j!) phi(x) integral_-inf^inf f(z) He_j(z) dz.}
#' The cumulative distribution function is 'simply' the integral of this
#' expansion. Abusing certain facts regarding the PDF and CDF of the normal
#' distribution and the probabilist's Hermite polynomials, the CDF has
#' the representation
#' \deqn{F(x) = \Phi(x) - \sum_{1\le j} \frac{He_{j-1}(x)}{j!} \phi(x) \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{F(x) = Phi(x) - sum_{1 <= j} (He_{j-1}(x)/j!) phi(x) integral_-inf^inf f(z) He_j(z) dz.}
#'
#' These series contain coefficients defined by the probability distribution 
#' under consideration. They take the form
#' \deqn{c_j \int_{-\infty}^{\infty} f(z) He_j(z) \mathrm{d}z.}{c_j = integral_-inf^inf f(z) He_j(z) dz.}
#' Using linearity of the integral, these coefficients are easily computed in
#' terms of the coefficients of the Hermite polynomials and the raw, uncentered
#' moments of the probability distribution under consideration. It may be the
#' case that the computation of these coefficients suffers from bad numerical
#' cancellation for some distributions, and that an alternative formulation
#' may be more numerically robust.
#'
#' @usage
#'
#' dapx.gca(x, raw.moments)
#'
#' papx.gca(q, raw.moments, lower.tail=TRUE)
#'
#' @param x where to evaluate the approximate density.
#' @param q where to evaluate the approximate distribution.
#' @param raw.moments an atomic array of the zeroth through kth raw moments
#' of the probability distribution. The first element should be a 1.
#' 2FIX: that seems like a stupid, arbitrary choice.
#' @param lower.tail whether to compute the lower tail. If false, we approximate the survival function.
#' @return The approximate density at \code{x}.
#'
#' @keywords distribution 
#' @seealso \code{\link{qapx.cf}}
#' @export 
#' @template ref-Jaschke
#' @aliases papx.gca 
#' @references
#'
#' S. Blinnikov and R. Moessner. "Expansions for nearly Gaussian
#' distributions." Astronomy and Astrophysics Supplement 130 (1998): 193-205.
#' \url{http://arxiv.org/abs/astro-ph/9711239}
#'
#' @examples 
#' # normal distribution:
#' xvals <- seq(-2,2,length.out=501)
#' d1 <- dapx.gca(xvals, c(1,0,1,0,3,0))
#' d2 <- dnorm(xvals)
#' d1 - d2
#'
#' qvals <- seq(-2,2,length.out=501)
#' p1 <- papx.gca(qvals, c(1,0,1,0,3,0))
#' p2 <- pnorm(qvals)
#' p1 - p2
#' @template etc
dapx.gca <- function(x,raw.moments) {
	raw.moments[1] <- 1  # just in case
	order.max <- length(raw.moments) - 1

	mu.central <- moments::raw2central(raw.moments)
	mu.std <- central2std(mu.central)
	eta <- (x - raw.moments[2]) / sqrt(mu.central[3])
	hermi <- orthopolynom::hermite.he.polynomials(order.max+1, normalized=FALSE)

	retval <- dnorm(eta)
	phi.eta <- retval
	for (iii in c(1:order.max)) {
		ci <- (sum(coef(hermi[[iii+1]]) * mu.std[1:(iii+1)])) / factorial(iii)
		retval <- retval - ci * phi.eta * as.function(hermi[[iii+1]])(eta)
	}

	# adjust back from standardized
	retval <- retval / sqrt(mu.central[3])
	return(retval)
}
#' @export 
papx.gca <- function(q,raw.moments,lower.tail=TRUE) {
	order.max <- length(raw.moments) - 1
	if (!lower.tail) {
		# transform q and the raw moments
		q <- - q;
		raw.moments <- raw.moments * (-1^(0:order.max))
	}

	mu.central <- moments::raw2central(raw.moments)
	mu.std <- central2std(mu.central)
	eta <- (q - raw.moments[2]) / sqrt(mu.central[3])
	hermi <- orthopolynom::hermite.he.polynomials(order.max, normalized=FALSE)

	retval <- pnorm(eta)
	phi.eta <- dnorm(eta)
	for (iii in c(1:order.max)) {
		ci <- (sum(coef(hermi[[iii+1]]) * mu.std[1:(iii+1)])) / factorial(iii)
		retval <- retval - ci * phi.eta * as.function(hermi[[iii]])(eta)
	}
	return(retval)
}

#' @title Approximate quantile via Cornish-Fisher expansion.
#'
#' @description 
#'
#' Approximate the quantile function of a distribution via its cumulants.
#'
#' @details
#'
#' Given the cumulants of a probability distribution, we approximate the 
#' quantile function via a Cornish-Fisher expansion.
#'
#' @usage
#'
#' qapx.cf(p, raw.cumulants)
#'
#' @param p where to evaluate the approximate distribution.
#' @param raw.cumulants an atomic array of the zeroth through kth raw cumulants. The first 
#' value is ignored, the second should be the mean of the distribution, the third should
#' be the variance of the distribution, the remainder are raw cumulants.
#' @return The approximate quantile at \code{x}.
#'
#' @keywords distribution 
#' @seealso \code{\link{dapx.gca}, \link{papx.gca}, \link{AS269}}
#' @export 
#' @template ref-AS269
#' @template ref-Jaschke
#' @examples 
#' # normal distribution:
#' pvals <- seq(0.001,0.999,length.out=501)
#' q1 <- qapx.cf(pvals, c(1,0,1,0,0,0,0,0))
#' q2 <- qnorm(pvals)
#' q1 - q2
#' @template etc
qapx.cf <- function(p,raw.cumulants) {
	order.max <- length(raw.cumulants) - 1
	# this should be a standard routine:
	std.cumulants <- raw.cumulants / (raw.cumulants[3] ^ ((0:(length(raw.cumulants)-1))/2))
	if (length(std.cumulants) > 3) {
		gammas <- std.cumulants[4:length(std.cumulants)]
		w <- AS269(z=qnorm(p),cumul=gammas,all.ords=FALSE)
	} else {
		w <- qnorm(p)
	}

	# now convert back from standardized:
	retval <- raw.cumulants[2] + w * sqrt(raw.cumulants[3])
	return(retval)
}



#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
