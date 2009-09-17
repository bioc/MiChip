\name{setIntensityCutoff}
\alias{setIntensityCutoff}
\title{Sets a cutoff for the lowest intensity value}
\description{
	Any value less than the cutoff value will be set to NA. This allows near background
	intensity values to be excluded
}
\usage{
	setIntensityCutoff(dmat, intensityCutoff)
}
\arguments{
	\item{dmat}{matrix of intensity values to which the cutoff value is applied}
	\item{intensityCutoff}{value of lowest acceptable intensity value in the experiment}
}
\examples{
#Set all the values under 50 in a matrix to NA
\dontrun{
dmatOver50 <- setIntensityCutoff(dmat, 50)
}
}
\keyword{data}
