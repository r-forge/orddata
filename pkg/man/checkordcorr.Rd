\name{check.ordcorr}
\alias{check.ordcorr}

\title{
Check Correlations for Ordinal Values
}
\description{
The function checks if the correlations matrix and the corresponding
probabilities lead to a result for \code{rmvord}.
}
\usage{
check.ordcorr(probs, Cor, switchmax = 1e+06, n = 100, showcommonprob = FALSE, showX = FALSE)
}

\arguments{
  \item{probs}{
    List of probabilities for the variables, length of probability equals
  number of items, length of list equals number of variables
}
  \item{Cor}{
Correlation matrix
}
  \item{switchmax}{
Number of maximal switches
}
  \item{n}{
Number of Observations
}
  \item{showcommonprob}{
If TRUE common probabilities are shown
}
  \item{showX}{
If TRUE resulting correlation matrix is shown
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a \code{logical} indicating if combination of \code{n}, \code{prob} and \code{Cor} is possible.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Dominik Traeger and Sebastian Kaiser
}

\seealso{
\code{\link{minmax.ordcorr}}
}
\examples{
check.ordcorr(list(c(1,1,1,1)/4,c(1,1,1,1)/4), cbind(c(0.5, 0.4), c(0.4, 0.8)), switchmax = 10000, n = 100, showcommonprob = TRUE, showX = TRUE)
}

\keyword{distribution}
\keyword{multivariate}