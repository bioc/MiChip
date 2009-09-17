\name{myForgivingMedian}
\alias{myForgivingMedian}
\title{Produce Median from Probe Intensity values}
\description{
	Creates a median to summarize the intensities for individual probes, giving that not
	all probes will have a valid intensity reading
}
\usage{
	myForgivingMedian(mat, minSumlength=0)
}
\arguments{
	\item{mat}{matrix of data to calculate the median from}
	\item{minSumlength}{The lowest acceptable length of the matrix to calculate a median}
}
\examples{

#Calculate the median of a matrix omiting NAs
\dontrun{
myForgivingMedian(mat, minSumlength=0)
}
}
\keyword{datagen}
