\name{boxplotDataNoFile}
\alias{boxplotDataNoFile}
\title{Create Boxplot of data}
\description{
	Creates a boxplot of expression data contained in a matrix.
}
\usage{
boxplotDataNoFile(dmat, exptname, dlevel)
}
\arguments{
	\item{dmat}{matrix containing expression data to be boxplotted}
	\item{exptname}{Name of the experiment, used to build plot title}
	\item{dlevel}{Stage of the experiment e.g. raw, summarized, normalized}
}
\examples{
#Create a boxplot containing normalized expression data for myexpt
\dontrun{boxplotDataNoFile(dmat, "MyExpt", "mednormed")}
}
\keyword{aplot}
