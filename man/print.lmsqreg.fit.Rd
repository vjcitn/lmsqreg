\name{print.lmsqreg.fit}
\alias{print.lmsqreg.fit}
\title{
Print information on a quantile regression
}
\description{
Print information on a quantile regression
}
\usage{
print.lmsqreg.fit(x,...)
}
\arguments{
\item{x}{
An object of class "lmsqreg.fit"
}}
\section{Side Effects}{
Prints information on convergence and
quality of an lmsqreg.fit.
}
\details{
Nominal and actual percentiles are compared; the
latter are computed by linear interpolation
between successive xtarget values across the range of x.
The adequacy of the transformations is also examined
within ranges of the independent variable: the Z scores
from the fitted model are subjected to the Kolmogorov Smirnov
test for standard normality.
}
\examples{


#lms quantile regression, version 2.4, fit date Tue Nov 12 19:25:55 EST 1996
#
#Dependent variable: jjj , independent variable: nnn 
#The fit converged with EDF=( 3,5,3 ), PL= 299.045 
#                                                             
#  nominal percentile 0.050 0.100 0.25 0.500 0.75 0.900 0.950
#estimated percentile 0.053 0.113 0.25 0.503 0.73 0.887 0.947
# 
#Shapiro Wilk tests: (intervals in nnn //p-values)
#  9.999+ thru 11.974 11.974+ thru 14.257 14.257+ thru 16.251 
#               0.239               0.334               0.568
# 16.251+ thru 17.915 17.915+ thru 19.965 Overall 
#               0.174               0.191   0.007
# 
#t tests: (intervals in nnn //p-values)
#  9.999+ thru 11.974 11.974+ thru 14.257 14.257+ thru 16.251 
#               0.185               0.025               0.548
# 16.251+ thru 17.915 17.915+ thru 19.965 Overall 
#               0.161               0.461   0.883
# 
#X2 tests (unit variance): (intervals in nnn //p-values)
#  9.999+ thru 11.974 11.974+ thru 14.257 14.257+ thru 16.251 
#               0.475               0.158               0.922
# 16.251+ thru 17.915 17.915+ thru 19.965 Overall 
#               0.266               0.808   0.809
#

}
\keyword{models}
\section{Version:}{
Document version 2.4 97/03/26
/usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.print.lmsqreg.fit.d
}
% Converted by Sd2Rd version 1.21.


