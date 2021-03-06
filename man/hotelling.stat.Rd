% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hotelling.test.r
\name{hotelling.stat}
\alias{hotelling.stat}
\alias{hotel.stat}
\title{Calculate Hotelling's two sample T-squared test statistic}
\usage{
hotelling.stat(x, y, shrinkage = FALSE)
}
\arguments{
\item{x}{a nx by p matrix containing the data points from sample 1}

\item{y}{a ny by p matrix containg the data points from sample 2}

\item{shrinkage}{set to \code{TRUE} if the covariance matrices are to be estimated
using Schaefer and Strimmer's James-Stein shrinkage estimator}
}
\value{
A list containing the following components:
\item{statistic}{Hotelling's (unscaled) T-squared statistic} \item{m}{The
scaling factor - this can be used by by multiplying it with the test
statistic, or dividing the critical F value} \item{df}{a vector of length
containing the numerator and denominator degrees of freedom} \item{nx}{The
sample size of sample 1} \item{ny}{The sample size of sample 2} \item{p}{The
number of variables to be used in the comparison}
}
\description{
Calculate Hotelling's T-squared test statistic for the difference in two
multivariate means.
}
\details{
Note, the sample size requirements are that nx + ny - 1 > p. The procedure
will stop if this is not met and the shrinkage estimator is not being used.
The shrinkage estimator has not been rigorously tested for this application
(small p, smaller n).
}
\examples{

data(container.df)
split.data = split(container.df[,-1],container.df$gp)
x = split.data[[1]]
y = split.data[[2]]
hotelling.stat(x, y)
hotelling.stat(x, y, TRUE)

}
\references{
Hotelling, H. (1931). ``The generalization of Student's ratio.''
Annals of Mathematical Statistics 2 (3): 360--378.

Schaefer, J., and K. Strimmer (2005). ``A shrinkage approach to large-scale
covariance matrix estimation and implications for functional genomics.''
Statist. Appl. Genet. Mol. Biol. 4: 32.

Opgen-Rhein, R., and K. Strimmer (2007). ``Accurate ranking of
differentially expressed genes by a distribution-free shrinkage approach.''
Statist. Appl. Genet. Mol. Biol. 6: 9.
}
\author{
James M. Curran
}
\keyword{htest}
