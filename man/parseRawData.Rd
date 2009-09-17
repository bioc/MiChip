\name{parseRawData}
\alias{parseRawData}
\title{Parse raw data files to create an ExpressionSet}
\description{
	Loads all the gpr scanner output files in a particular directory and returns an 
	\code{ExpressionSet} of the hybridizations in a MiChip experiment
}
\usage{
parseRawData(datadir=".", pat="gpr")
}
\arguments{
	\item{datadir}{a directory containing one or my files of scanner output 
	from MiChip hybridizations}
	\item{pat}{a string containing the three letter extension of the scanner output files}
}

\examples{
\dontrun{
## Load all *.gpr files in current directory
parseRawData(datadir=".", pat="gpr")


## Load all *.gpr files in a specified directory , windows
parseRawData(datadir="c:\\mydata\\grpdata\\expt1\\", pat="gpr")


## Load all *.gpr files in a specified directory, linux
parseRawData(datadir="/home/myuser/gprdata/extp1/", pat="gpr")
}
}
\keyword{file}
