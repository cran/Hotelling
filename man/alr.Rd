\name{alr}
\alias{alr}
\title{Additive log ratio transformation}
\description{Aitchison's additive log ratio tranformation for
  compositional data
}
\usage{
alr(form, data, group = NULL)
}
\arguments{
  \item{form}{a formula which specifies the denominator variable as the
    response}
  \item{data}{a data frame in which the data is stored}
  \item{group}{if not NULL then a character string specifying the name
    of the grouping variable}
}
\value{
  a data frame with the ALR transformation applied to data. Each row in
  the data frame is standardized with respect to a specific variable by
  dividing by that variable. The logarithms of the resulting ratios are returned. If a grouping variable is specified, then this is preserved.
}
\details{
  This function will give a warning if zeros are present because the
  transformed data will have -Infs.
}
\references{
  Aitchison, J. (1986). ``The Statistical Analysis of Compositional Data'', Chapman and Hall, reprinted in 2003 with additional material by The Blackburn Press
}
\seealso{
  clr
}
\examples{
data(bottle.df)

## transform with respect to manganese
alr(Mn~., bottle.df, "Number")

## transform the data with respect to barium, but removing the
## grouping in column 1
alr(Ba~., bottle.df[,-1])
}
\author{James M. Curran}
\keyword{transformation}
