\name{validate.report}
\alias{validate.report}
\title{
LMS quantile regression stats on test data
}
\description{
report on fit of quantile regression model to arbitrary data
}
\usage{
validate.report(qreg, y.val, t.val, rule=2, xcuts=quantile(t.val, c(0.2, 0.4, 0.6, 0.8)))
}
\arguments{
\item{qreg}{
an lmsqreg.fit object
}
\item{y.val}{
values of dependent variable
}
\item{t.val}{
values of independent variable
}
\item{rule}{
interpolation rule parameter
}
\item{xcuts}{
levels of independent variable within which
goodness of fit should be assessed
}}
\value{
prints out goodness of fit appraisal for data (y,t)
which were presumably not used in the fit that produced qreg.
}
\keyword{models}
% Converted by Sd2Rd version 1.21.


