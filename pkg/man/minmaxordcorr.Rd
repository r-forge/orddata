\name{minmax.ordcorr}
\alias{minmax.ordcorr}

\title{
Minimum and Maximum of Possible Correlations
}
\description{
  This program \code{minmax.ordcorr} checks the first condition of the
  feasibility of a correlation matrix of ordinal random numbers. A mean vector (as list) needs to be specified. It returns yes/no if also a correlation matrix
  was given and in either case the Min-Max Correlation Matrix, which has
  the minimum correlation in the lower triangular matrix and the maximum
  correlation in the upper triangular matrix.
}
\usage{
minmax.ordcorr(probs, Cor = 0, n = 1e+06, showX = FALSE)
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
  \item{showX}{
If TRUE resulting correlation matrix is shown
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  It returns yes/no if also a correlation matrix was given and
  in either case the Min-Max Correlation Matrix, which has
  the minimum correlation in the lower triangular matrix and the maximum
  correlation in the upper triangular matrix.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Dominik Traeger and Sebastian Kaiser
}

\seealso{
\code{\link{check.ordcorr}}
}
\examples{
minmax.ordcorr(list(c(1,1,1,1)/4,c(1,1,1,1)/4), cbind(c(0.5, 0.4), c(0.4, 0.8)), n = 1000, showX = TRUE)
}

\keyword{distribution}
\keyword{multivariate}