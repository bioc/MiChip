\name{returnAnnotatedDataMatrix}
\alias{returnAnnotatedDataMatrix}
\title{returns and annotated data matrix from an ExpressionSet}
\description{
	Takes an \code{ExpressionSet} and returns a data martix of 
	feature annotation to the left and hyb specific expression/flag data to the right
}
\usage{
returnAnnotatedDataMatrix(eset, dataElement)
}
\arguments{
	\item{eset}{\code{ExpressionSet} containing the matrix and annotation to output}
	\item{dataElement}{a string containing the name of the data element in the ExpressionSet
		to be output}
}
\examples{

#Write out an annotated tab delimited file for the normalized data
\dontrun{
returnAnnotatedDataMatrix(normedEset,"exprs")}
}
\keyword{file}
