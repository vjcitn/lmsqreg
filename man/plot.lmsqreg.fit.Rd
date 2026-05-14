\name{plot.lmsqreg.fit}
\alias{plot.lmsqreg.fit}
\title{
Plot a quantile regression
}
\description{
Plot centiles and transformation functions
}
\usage{
plot.lmsqreg.fit(x, fullsys=T, medname="P0.5", xlab=" ", ylab=" ", Yshift=0, 
     title=paste("LMS fit with edf = 
     (", ob$edf[1], ",", ob$edf[2], ",", ob$edf[3], ")", 
     CEX = 1.1, tx = function(z) z,
    xaxat = NULL, xaxlab = NULL, ...)

}
\arguments{
\item{obj}{
An object of class "lmsqreg.fit"
}
\item{fullsys}{
If T, plot quantile system and L, M, S functions on 
four currently valid screens; otherwise
just plot quantile system on the current valid screen.
The function will run split.screen if there are
not exactly four currently valid screens.
}
\item{medname}{
The character variable naming the row where
the median estimates are kept in obj[["qsys"]].
}
\item{xlab}{
Label for plot x-axis.
}
\item{ylab}{
Label for plot y-axis.
}
\item{title}{
Title of plot
}
\item{tx}{function to transform x}
\item{xaxat}{points at which ticks should be drawn on x axis}
\item{xaxlab}{labels at ticks drawn on x axis}
}
\details{
Currently uses layout to set up a display of four graphs.
}
\section{Side Effects}{
Uses layout.
}
\keyword{models}
\section{Version info}{
Document version 2.4 97/03/26
/usr16/stdevs/stdev0f/SLIBS/lmsqreg.dev.obs/SCCS/s.plot.lmsqreg.fit.d
}
% Converted by Sd2Rd version 1.21.

