\name{rmvord_mc}
\alias{rmvord_mc}
\title{Multivariate Ordinal Random Variates by Monte Carlo Simulation}

\description{
Creates correlated multivariate ordinal random variables by a Monte
Carlo simulation.
}
\usage{
rmvord_mc(probs, Cor, n = 1)
}
%- maybe also 'usage' for other objects documented here.
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
}

\details{
Ordinal values are produced by shifting the variables until correlation
structure is reached.
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
\code{\link{rmvord_b}},\code{\link{rmvord}}
}
\examples{
rmvord_mc(list(c(1,1,1,1)/4,c(1,1,1,1)/4), cbind(c(1, 0.4), c(0.4, 1)), n = 20)
}

\keyword{distribution}
\keyword{multivariate}

