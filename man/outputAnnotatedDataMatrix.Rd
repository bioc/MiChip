\name{outputAnnotatedDataMatrix}
\alias{outputAnnotatedDataMatrix}
\title{Outputs a tab delimited file from an ExpressionSet}
\description{
	Takes an \code{ExpressionSet} and outputs a tab delimited file containing
	feature annotation to the left and hyb specific expression/flag data to the right
}
\usage{
outputAnnotatedDataMatrix(eset, exptname, stage, dataElement)
}
\arguments{
	\item{eset}{\code{ExpressionSet} containing the matrix and annotation to output}
	\item{exptname}{a string containing the name of the experiment. Used to build file name}
	\item{stage}{a string containing the stage of the data in the matrix e.g. normalized}
	\item{dataElement}{a string containing the name of the data element in the ExpressionSet
		to be output}
}
\examples{

#Write out an annotated tab delimited file for the normalized data
\dontrun{
outputAnnotatedDataMatrix(normedEset, "MyMicroArrayExpt","Median_Normalized", "exprs")}
}
\keyword{file}
