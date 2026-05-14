\name{lmsqreg.object}
\docType{data}
\alias{lmsqreg.object}
\title{ Object representing an LMS quantile regression }
\description{ Output of lmsqreg.fit, an object of class "lmsqreg.fit".  }

\section{GENERATION:}{ Such objects are produced by calls to lmsqreg.fit.  }
\section{METHODS:}{ Currently methods for printing and plotting these objects exist.  }
\section{lms.ans:}{
There are three main lists in any lmsqreg.fit object.
The first, lms.ans, is a list including components

\describe{
\item{ordt}{ the ordered abscissa values on which the fit is based }
\item{lam, mu, sig}{ the estimated lambda, mu, sigma function values corresponding to ordt }
\item{edf}{ requested equivalent degrees of freedom }
\item{niter}{ actual number of iterations used in fit }
\item{converged}{ logical variable indicating whether or not fit converged }
\item{fit.date}{ date this fit was created }
\item{fitter.version}{ version of the fitting software }
\item{yname}{ name of S object used as dependent variable }
\item{xname}{ name of S object used as independent variable The second component, qsys, is a list including components
}
\item{outmat}{ a matrix giving the estimated centile function values }
\item{targetx}{ a vector giving the abscissa values at which centile estimates were requested }
\item{pcts}{ vector giving requested centiles }
\item{edf}{ requested equivalent degrees of freedom }
}
}

\section{validout:}{
\describe{
\item{validout}{
The third component, validout, is a list with
two components: the first is a vector giving the
requested (nominal) centile values C, the second
is the vector of estimated percentile values --
the proportion of data lying beneath the Cth centile
function.
}
}
}

\note{Version: Document version 2.3 97/03/26 /usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.lmsqreg.object.d
}
\keyword{models}

% Converted by Sd2Rd version 1.21.


