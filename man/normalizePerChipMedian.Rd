\name{normalizePerChipMedian}
\alias{normalizePerChipMedian}
\title{Normalize to median intensity}
\description{
	Normalizes intensity values to the median of each chip
}
\usage{
normalizePerChipMedian(eset)

}
\arguments{
	\item{eset}{\code{ExpressionSet} containing chip intensity values to be normalized}
}
\examples{
#Normalize expression data in an Eset to the median
\dontrun{
normedDataEset <- normalizePerChipMedian(eset)
}
}
\keyword{models}
