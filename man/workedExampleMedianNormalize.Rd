\name{workedExampleMedianNormalize}
\alias{workedExampleMedianNormalize}
\title{Worked Example of MiChip Processing}
\description{
	Loads a set of hybridizations into a matrix and them proceeds to filter, summarize and 
	median normalize them
}
\usage{
workedExampleMedianNormalize(exptname, intensityCutoff=0, datadir=".", minSumlength, madAdjust = FALSE)
}
\arguments{
	\item{exptname}{string indicating the name of the experiment}
	\item{intensityCutoff}{The intensity value for accepting the spots intensity value in the experiment}
	\item{datadir}{The directory where hybridization files are found.}
	\item{minSumlength}{Minimum exceptable number of values to summarize intensity value.}
	\item{madAdjust}{if \code{TRUE} then summarized data will be filtered according to the MAD median absolute deviation and set to  \code{NA} if the median is less than MAD}
}
\examples{
#Normalize data in the current directory to the median per chip
datadir <- system.file("extdata", package="MiChip")
myNormedEset <-workedExampleMedianNormalize("MyExpt", intensityCutoff=0, datadir, minSumlength=0, madAdjust=TRUE)
}
\keyword{design}
