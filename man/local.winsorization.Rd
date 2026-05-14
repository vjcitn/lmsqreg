\name{local.winsorization}
\alias{local.winsorization}
\alias{locwin.envelope}
\title{
Relocation of ``locally outlying'' data points 
}
\description{
For use with LMS quantile regression
}
\usage{
local.winsorization(x, y, ncut=5, k=20)
~move the above line to just above the first optional argument
}
\arguments{
\item{x}{
independent variable 
}
\item{y}{
dependent variable, same length as x
}
\item{ncut}{
number of equal-sized groups formed on values
of x within which outliers are sought (relative to
group mean, using GESD (Rosner, Technometrics 1983)
}
\item{k}{
upper bound on number of outliers suspected in any group
}}
\value{
a list with components x, y and bad -- where bad gives
the indices of any outliers, y[bad] are winsorized
values, x and y[-bad] are the original data values
in original order
}
\details{
local winsorization is preferable to outlier deletion
for quantile regression procedures; LMS is particularly
sensitive to aberrant values and can fail to converge
or even to iterate in the presence of outliers.
}
\references{
}
\examples{
     set.seed(123)
      nnn <- runif(300, 10, 20)
      jjj <- 8 + 2 * sin(nnn) + rnorm(100, 0, nnn/11)
      jjj[1] <- jjj[1]+100
      tfff <- try(fff <- lmsqreg.fit(jjj, nnn))
      if (!inherits(tfff,"try-error")) print(fff)
      LW <- local.winsorization(nnn,jjj)
      fxjjj <- LW$y
      fffx <- lmsqreg.fit(fxjjj, nnn)
      fffx
}

\keyword{outliers}
\keyword{model}
% Converted by Sd2Rd version 1.21.


