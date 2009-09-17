\name{workedExampleNotNormalizedData}
\alias{workedExampleNotNormalizedData}
\title{Worked Example of MiChip Processing}
\description{
	Loads a set of hybridizations into a matrix and them proceeds to filter and summarize these data
}
\usage{
workedExampleNotNormalizedData(exptname, intensityCutoff=0, datadir=".", minSumlength, madAdjust = FALSE )
}
\arguments{
	\item{exptname}{string indicating the name of the experiment}
	\item{intensityCutoff}{The intensity value for accepting the spots intensity value in the experment}
	\item{datadir}{The directory contain data from the experiment}
	\item{minSumlength}{Minimum exceptable number of values to summarize intensity value.}
	\item{madAdjust}{if \code{TRUE} then summarized data will be filtered according to the MAD median absolute deviation and set to  \code{NA} if the median is less than MAD}

}
\examples{
#Summarizes the data in the current directory 
\dontrun{
mySummarizedEset <-workedExampleNotNormalizedData("MyExpt", intensityCutoff=0, datadir=".", minSumlength=0, madAdjust=TRUE)
}
}
\keyword{design}
