\name{standardRemoveRows}
\alias{standardRemoveRows}
\title{Removes a standard list of features for MiChip processing}
\description{
	Removes all empty spots, control spots, U6 RNA, non human spots from an 
	\code{ExpressionSet} in the standard fashion. A wrapper for removeUnwantedRows
}
\usage{
	standardRemoveRows(rawData)
}
\arguments{
	\item{rawData}{\code{ExpressionSet} containing the matrix to be filtered}
}
\examples{
#Filter standard rows from an ExpressionSet
\dontrun{
myfilterdEset <-standardRemoveRows(rawData)
}
}
\keyword{data}
