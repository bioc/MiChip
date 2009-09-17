\name{removeUnwantedRows}
\alias{removeUnwantedRows}
\title{Removes unwanted rows from data matrix}
\description{
	Due to the requirements of spotting the chips, some of the
	spots are empty. Others contain controls or features from
	another species that may not be wanted in the analysis.
	This method removes them
}
\usage{
removeUnwantedRows(rawData, filters)
}
\arguments{
	\item{rawData}{\code{ExpressionSet} containing matrix of data to be filtered}
	\item{filters}{list of strings to be to be filtered from annotation gene name column}
}
\examples{
#Removes empty and control spots from data matrix
\dontrun{
filters=c("empty", "control")
filteredData <- removeUnwantedRows(rawData, filters)}
}
\keyword{data}
