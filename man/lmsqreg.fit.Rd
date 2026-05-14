\name{lmsqreg.fit}
\alias{lmsqreg.fit}
\alias{gam.spar.R}
\alias{get.quantiles}
\alias{lmsqreg.search}

\title{
Fit A Quantile Regression
}
\description{
Cole and Green LMS (Stat in Med, 1992)
}
\usage{
lmsqreg.fit(YY, TT, edf = c(3, 5, 3), targlen = 50, targetx = seq(min(TT), max(TT),
        length = targlen), pvec = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), 
        maxit = 15, tol = 0.01, verb = FALSE, lam.fixed = NULL, mu.fixed = NULL, 
        sig.fixed = NULL, xcuts = quantile(TT, c(0.2, 0.4, 0.6, 0.8)),
sig.init, lam.init)
}
\arguments{
\item{YY}{
Ordinate values; must be positive
}
\item{TT}{
Abscissa values
}
\item{edf}{
A 3-vector specifying
"equivalent degrees of freedom" for
Box-Cox transformation, median, and standard deviation
functions respectively, assumed to be smooth in TT.  The scale
of edf corresponds to the df parameter of s() in
gam(), for which df=1 corresponds to a linear model.
Constant models for the component
functions can be obtained
by setting lam.fixed, etc. to specified values (see below).
}
\item{targlen}{
Number of points at which smooth estimates of L, M, S
should be extracted for quantile plotting; if quantile
plots are jagged in appearance, the value of this parameter should
be increased.  This parameter has no effect on the
fitting process.
}
\item{targetx}{
Points on which smooth estimates of L, M, S should
be extracted for quantile plotting
}
\item{pvec}{
Vector of target percentiles for plotting; default
(5,10,25,50,75,90,95) percentiles will be plotted.
}
\item{maxit}{
Limit to number of Fisher scoring iterations
}
\item{tol}{
Tolerance on change between estimates of L, M, S
on successive scoring iterations
}
\item{verb}{
verbose run; will give iteration-specific information if T
}
\item{lam.fixed}{
if NULL, lambda will be estimated; if non-null (numeric atom),
lambda will be set to this value.  Set to zero to force
a log transformation; set to unity to perform no transformation.
}
\item{mu.fixed}{
if NULL, mu will be estimated; if non-null (numeric atom),
mu will be set to this value.
}
\item{sig.fixed}{
if NULL, sigma (coefficient of variation function) 
will be estimated; if non-null (numeric atom),
sigma will be set to this value.
}
\item{xcuts}{
vector of x values defining classes within which Kolmogorov Smirnov (KS)
tests for normality of derived Z-scores will be conducted.
}
\item{sig.init}{initialization}
\item{lam.init}{initialization}
}
\value{
A list of class "lmsqreg.fit"; see lmsqreg.object for details.
}
\details{
}
\references{
TJ Cole, PJ Green,
Smoothing reference centile curves: The LMS method and penalized likelihood,
Statistics in Medicine, v11, 1992, p1305--1319.
}
\seealso{
lmsqreg.object, print.lmsqreg.fit, plot.lmsqreg.fit
}
\examples{
 set.seed(123)
 nnn <- runif(300, 10, 20)
 jjj <- 8 + 2 * sin(nnn) + rnorm(100, 0, nnn/11)
 fff <- lmsqreg.fit(jjj, nnn)
 fff      
# results from SPLUS:
#lms quantile regression, version 2.3, fit date Sat May 25 20:44:57 EDT 1996
#
# 
#Dependent variable: jjj , independent variable: nnn 
#The fit converged with EDF=( 3,5,3 ), PL= 299.045 
#
#                                                             
#  nominal percentile 0.050 0.100 0.25 0.500 0.75 0.900 0.950
#estimated percentile 0.053 0.113 0.25 0.503 0.73 0.887 0.947
#
# 
#Shapiro Wilk tests: (intervals in nnn //p-values)
#  9.999+ thru 11.974 11.974+ thru 14.257 14.257+ thru 16.251 
#               0.239               0.334               0.568
# 16.251+ thru 17.915 17.915+ thru 19.965 
#               0.174               0.191
#

}
\keyword{quantile}
\keyword{models}
\section{Version information:}{
Document version 2.4 97/03/26
/usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.lmsqreg.fit.d
}
% Converted by Sd2Rd version 1.21.


