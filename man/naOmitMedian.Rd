\name{naOmitMedian}
\alias{naOmitMedian}
\title{Calculates the median of an array excluding NAs}
\description{
	Calculates the median of an array excluding NAs
	
}
\usage{
naOmitMedian(mat, madAdjust=FALSE)
}
\arguments{
	\item{mat}{A single dimensional matrix}
	\item{madAdjust}{if \code{TRUE} then summarized data will be filtered according to the MAD median absolute deviation and set to  \code{NA} if the median is less than MAD}

}
\examples{
	#Calculate the median of matrix mat omiting NAs
\dontrun{
	myMedian <-naOmitMedian(mat, madAdjust=TRUE)
}
}
\keyword{datagen}
