\name{zscores}
\alias{zscores}
\title{
zscores from an lmsqreg system
}
\description{
use quantile regression object to generate zscores for new data
}
\usage{
zscores(y, x, obj)
}
\arguments{
\item{y}{
y values for which z-scores are desired
}
\item{x}{
x values corresponding to y values
}
\item{obj}{
an object of class lmsqreg.fit
}}
\value{
A vector of z-scores (signed number of standard deviations between y and 
its median value given the corresponding x
as determined by the lmsqreg.fit object.  Specifically,
the lmsqreg.fit object determines functions M(x) = med(Y|X=x)
and S(x) = S.D.(Y|X=x), and the zscore for the observation (y,x)
is (y-M(x))/S(x).
}
\details{
To construct M(x) and S(x),
linear interpolation is used when x lies between x values
on which obj was fit; constant extrapolation is used when
x lies out of the range of x used to create obj.  In the
latter case, warning messages are produced with information
on the magnitude of the discrepancy between x and the
range of x on which obj was based.
}
\references{
TJ Cole, PJ Green, Smoothing reference centile curves: The
LMS method and penalized likelihood; Stat in Med, 1992 (p1305).
}
\seealso{
lmsqreg.fit
}
\examples{
#> zscores(6,10,fff)
#[1] -0.03834572
#Warning messages:
#  constant extrap. from 10.0002 to 10 in: zscores(6, 10, fff)
## see lmsqreg.fit help for creation of fff; this example
## illustrates a minor extrapolation problem.
}
\keyword{quantile}
\keyword{models}
\section{Version information:}{
Document version 1.3 97/03/26
/usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.zscores.d
}
% Converted by Sd2Rd version 1.21.



