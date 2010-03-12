\name{rmvord_naiv}
\alias{rmvord_naiv}

\title{
Multivariate Ordinal Random Variates}

\description{
Creates correlated multivariate ordinal random variables by
  thresholding a normal distribution.
}
\usage{
rmvord_naiv(probs, Cor, n = 1, showCor = TRUE)
}

\arguments{
  \item{probs}{
    List of probabilities for the variables, length of probability equals
  number of items, length of list equals number of variables
}
  \item{Cor}{
Correlation matrix
}
  \item{n}{
Number of Observations
}
  \item{showCor}{
If TRUE correlation matrix is printed.
}

}
\details{
An analytical method was used to transform the correlation matrix.
}
\value{
Returns \code{n} observations with given marginal probabilities
\code{probs} and correlation structure \code{Cor}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Dominik Traeger and Sebastian Kaiser
}

\seealso{
\code{\link{rmvord_b}},\code{\link{rmvord_mc}}
}
\examples{
rmvord_naiv(list(c(1,1,1,1)/4,c(1,1,1,1)/4), cbind(c(1, 0.4), c(0.4, 1)), n = 20)
}

\keyword{distribution}
\keyword{multivariate}