\name{summarizeIntensitiesAsMedian }
\alias{summarizeIntensitiesAsMedian}
\title{Summarizes the probe intensity as median of replicates spotted}
\description{
	As the probes are spotted onto the in quaduplet
	or duplicate the values have to be combined in some way. This function
	takes the median of the intensities for the spots. Effectively
	the mean for duplicates. If less than half of the spots are
	present an NA is added
}
\usage{
	summarizeIntensitiesAsMedian(eset,minSumlength=0, madAdjust=FALSE)
}
\arguments{
	\item{eset}{ExpressionSet containing probe intensity data to be summarized}
	\item{minSumlength}{The lowest acceptable length of the matrix to calculate a median}
	\item{madAdjust}{if \code{TRUE} then summarized data will be filtered according to the MAD median absolute deviation and set to  \code{NA} if the median is less than MAD}
}
\examples{

#Calculate the median of a matrix omiting NAs
\dontrun{
summarizeIntensitiesAsMedian(eset,minSumlength=0,madAdjust=TRUE)
}
}
\keyword{datagen}
