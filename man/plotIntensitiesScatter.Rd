\name{plotIntensitiesScatter}
\alias{plotIntensitiesScatter}
\title{Plot pairwise intensity scatter}
\description{
	Creates a pairwise set of scatter plots from a data matrix and writes it out to file 
}
\usage{
plotIntensitiesScatter(dmat, controls=NULL, exptname, maintitle)
}
\arguments{
	\item{dmat}{matrix containing data from an experiment to be plotted}
	\item{controls}{matrix of row numbers containing control data to be plotted in a different colour}
	\item{exptname}{Name of the experiment, used for build the filename}
	\item{maintitle}{String used to build the maintitle of the graph}
}
\examples{
#Plot the pairwise intensities from myexpt
\dontrun{
plotIntensitiesScatter(dmat, NULL, "MyExpt", "Median_Normalized")
}
}
\keyword{hplot}
