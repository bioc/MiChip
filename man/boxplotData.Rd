\name{boxplotData}
\alias{boxplotData}
\title{Create Boxplot of data}
\description{
	Creates a boxplot of expression data contained in a matrix and writes this to a file.
}
\usage{
	boxplotData(dmat, exptname, dlevel)
}
\arguments{
	\item{dmat}{matrix containing expression data to be boxplotted}
	\item{exptname}{Name of the experiment, used to build filename}
	\item{dlevel}{Stage of the experiment e.g. raw, summarized, normalized}
}
\examples{

#Create a file of a boxplot containing normalized expression data for myexpt
\dontrun{
boxplotData(dmat, "MyExpt", "mednormed")
}
}
\keyword{aplot}

