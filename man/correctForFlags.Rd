\name{correctForFlags}
\alias{correctForFlags}
\title{Corrects for spots flagged as not present}
\description{
	Spots flagged with a -ve quality flag value by the scanner may be regarded as not present.
	This method sets their intensity to NA.
}
\usage{
correctForFlags(eset, intensityCutoff=0)
}
\arguments{
	\item{eset}{\code{ExpressionSet} containing intensity values and flags to be filtered}
	\item{intensityCutoff}{value of lowest acceptable intensity value in the experiment}	
}
\examples{
#Correct ExpressionSet for flags of spots marked as unreadable
\dontrun{
myCorrectedEset <-correctForFlags(eset, intensityCutoff=0)
}
}
\keyword{data}
