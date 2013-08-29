\name{print.hotelling.test}
\alias{print.hotelling.test}
\title{Prints the results from a Hotelling's
T-squared test for the difference in two multivariate sample means}
\description{
  Prints the test stastic, degrees of freedom and P-value
  from Hotelling's T-squared test for the difference in two multivariate
  sample means
}
\usage{
\method{print}{hotelling.test}(x,\dots)
}
\arguments{
  \item{x}{an object of type hotelling.test}
  \item{\dots}{any additional arguments to be passed to the hist command}
}
\examples{
data(bottle.df)
bottle.df = subset(bottle.df, Number == 1)
bottle.df$Number = rep(1:2,c(10,10))
fit = hotelling.test(.~Number, bottle.df, perm = TRUE)
fit
fit = hotelling.test(.~Number, bottle.df)
fit

## an explict call
print(fit)
}
\author{James M. Curran}
\keyword{print}
