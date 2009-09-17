\name{panelCor}
\alias{panelCor}
\title{Add Pearson Correlation value to plots}
\description{
	Adds a pearson correlation value to the scatter plots
}
\usage{
	panelCor(x,y, digits=2, prefix="r=")
}
\arguments{
	\item{x}{matrix of x values}
	\item{y}{matrix of y values to correlate with x}
	\item{digits}{number of digits to display}
	\item{prefix}{The string prefix that should be display on the scatterplot panel}
}
\examples{

#Calculate the median of a matrix omiting NAs
\dontrun{
panelCor(x,y, digits=2, prefix="r=")
}
}
\keyword{aplot}
