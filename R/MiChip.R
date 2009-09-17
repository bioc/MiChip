#############################################################
# ESetParseMichip.R
# Parses Out a set of genepix files to produce an ExpressionSet
# Date: 27.11.2007
# Author: Jonathon blake
# Version: 0.5.1
#############################################################
#Modifications
###############
#27.03.08 This version was used for incorporating in to MiChip
#library 0.02 and testing
#03.06.08 Version 0.4 version built removing the latest SNORD
#features from the chips.
#17.03.2009 Version 0.5 built to deal with changes in the 
#specification of names of the non representative probes and 
#the "" empties
#03.06.2009 Version 0.5.0 built dealing with -1 probe names
#and correcting grep errors
#29.06.2009 Version 0.5.1 Dealing with Bioconductor technical 
#review issues 
#02.07.2009 Version 0.5.2 Dealing with Bioconductor technical 
#review issues 2
#08.09.2009 Version 0.5.3 Dealing with Bioconductor peer review
#issues
#############################################################
library("Biobase")
#############################################################
#FUNCTIONS
#############################################################
#############################################################
#FILEHANDLING FUNCTIONS
#############################################################
################
#parseRawData
#Loads the data in a set of text files into an ExpressionSet for further use
#This method returns an ExpressionSet of median intensity - background and flag with for each
#spot and chip, The individual spots are identifiable by block row column.
#The ID column contains the probe ID which is spotted in quadruplicate (experimental)
#and possibly more often for controls
#The Name column contains the name for the microRNA which is supposedly represented by the probe
################
parseRawData<-function(datadir=".", pat="gpr")
{
	#Get the files in the directory specified
	chipfiles <- list.files(path=datadir, pattern=pat)
	pathfiles <- chipfiles
	
	if(datadir !=".") 
	{
		#pathfiles <-paste(datadir, pathfiles, sep="") 0.5.1
		pathfiles <- file.path(datadir, pathfiles) 
	
		#cat(pathfiles,"\n")
	}
	#Being verbose at the moment so check number of files
	if(length(chipfiles) == 0)
	{
		#Stop here there are no files to process
		stop("No files found matching file extension")
		
	}
	else
	{
		cat(length(chipfiles), "files found\n")
	
	}
	#Get the line number for the header line in these files
	
	#headlineno <- grep("Block", scan(chipfiles[1], nlines=40, quiet=TRUE, what="", sep="\n"))
#0.5.1	cat("Dealing with file ", pathfiles[1], "\n")
	headlineno <- grep("Block", scan(pathfiles[1], nlines=60, quiet=TRUE, what="", sep="\n"))
	if(!any(headlineno))
	{
		stop("No header found")
	}
	#Get the headerline
	#headercols <- scan(chipfiles[1], skip=headlineno-1, nlines=1, sep="\t", quiet=TRUE, what="")
#0.5.1	cat("Dealing with file ", pathfiles[1], "\n")
	headercols <- scan(pathfiles[1], skip=headlineno-1, nlines=1, sep="\t", quiet=TRUE, what="")
	#Get the column for the Name
	namecol <- grep("Name", headercols)
	#Get the column for the ID
	idcol <- grep("ID", headercols)
	#Get the column number for the Median - background
	medMinusBackcol <- grep("F532 Median - B532", headercols)
	#Get the column for flags
	flagscol <- grep("Flags", headercols)
	#Does anything have to be chomped?
	
	#Calculate the number of spots on the array
	ngenes <- length(scan(pathfiles[1], skip=headlineno, sep="\n", quiet=TRUE, what=""))
	#Read the first file to get the names of the features. The files must all be of the same
	#length otherwise this indicates different versions of the chip. 
	#DIFFERENT CHIP VERSIONS MUST THROW AN ERROR
	data <- read.table(pathfiles[1], skip=headlineno, sep="\t", header=FALSE)
	cat("Dealing with file ", chipfiles[1], "\n")
	coords<- c(paste(data[,1], "_", data[,2],"_", data[,3]))
	fd <- cbind(coords, as.data.frame(data[,idcol]), data[,namecol])
	colnames(fd) <- cbind("BlockID", "PROBEID", "NAME")
	feats <- new("AnnotatedDataFrame", as.data.frame(fd))
	#####
	#Build the meta data for the varlabels
	#####
	labdata <- as.data.frame(c("Composite of Spot Coordinates", "Exiqon Probe ID", "miRNA Name"))
	colnames(labdata) <- 'labelDescription'
	varMetadata(feats)<- labdata
	
	#####
	#Prepare matrices for the build
	#####
	#Extract the Sample Header info
	samples<- substring(chipfiles[1], 1, nchar(chipfiles[1])-4)
	#The flags data
	flags <- data[,flagscol]
	#The experimental intensity - background
	expdata <- data[,medMinusBackcol]
	
	#####
	#Go through the other files in the directory
	#####
	for(x in 2:length(chipfiles))
	{
		ngenesX <- length(scan(pathfiles[x], skip=headlineno, sep="\n", quiet=TRUE, what=""))
		if(ngenesX != ngenes)
		{
				stop("Different Number of Probes, looks like different versions of MiChip!")
		}
		data <- read.table(pathfiles[x], skip=headlineno, sep="\t", header=FALSE, quote="")
		
		cat("Dealing with file ", chipfiles[x], "\n")
		samples <- c(samples, substring(chipfiles[x], 1, nchar(chipfiles[x])-4))
		flags   <- cbind(flags,data[,flagscol])
		expdata <- cbind(expdata,data[,medMinusBackcol])
		
			
	}
	#Create an annotated dataframe for phenotype data, the sample names here
	#Add the column names to the dataframes and matrix
	sframe = as.data.frame(samples, stringsAsFactors=FALSE)
	colnames(sframe) <- "Samples"
	colnames(flags) = colnames(expdata) = rownames(sframe) = samples
	
	phenD <- new("AnnotatedDataFrame", data= sframe)
	
	#Build the expressionset
	eset <- new("ExpressionSet", exprs=expdata, featureData=feats, phenoData=phenD)
	#Add the flags to eset
	assayDataElement(eset, "flags") <- flags
	#Now have everything build return this as the  container of the raw data
	return(eset)
	
}
#############################################################
#############################################################
#Name: outputAnnotatedDataMatrix
#Description:Takes the expressionSet, and the name of an
#assayDataElement and outputs this into a text file with 
#the annotation provided in the featureData
#Author: Jonathon Blake
#Date:30.11.2007 
#############################################################
outputAnnotatedDataMatrix <-function(eset, exptname, stage, dataElement)
{
	filename <-  paste(exptname, "_",stage, ".txt", sep="")
	outtab <- returnAnnotatedDataMatrix(eset, dataElement)
	colarr <- colnames(outtab)
	write.table(as.data.frame(outtab), file=filename, sep="\t", quote=FALSE, col.names=colarr, row.names=FALSE)
	cat("File saved to ", file.path(getwd(), filename),"\n", sep="")
}
#############################################################
#############################################################
#Name: returnAnnotatedDataMatrix
#Description:Takes the expressionSet and returns an annotated
#data matrix
#Author: Jonathon Blake
#Date:08.09.2009
#############################################################
returnAnnotatedDataMatrix  <-function(eset, dataElement)
{
	dmat <- assayDataElement(eset, dataElement)
	pmat <- pData(featureData(eset))
	outtab <- cbind(pmat, dmat)
	return(outtab)
}
#############################################################
#DATA CLEANING FUNCTIONS
#############################################################
#############################################################
#Name:removeUnwantedRows
#Description: Cleans up the rawDataFrame removing unwanted rows
#such as the empty rows (approx half the chip) and unused controls
#We have a list of columns which should be removed normally, according
#to the spotting regime at EMBL. As such you would normally call the 
#standardRemoveRows function which is a wrapper for this function
#Call this method if you wish to remove a non-standard list of probe names.
#Author:Jonathon Blake
#Date:30.05.2007
#Date:28.11.2007: Updated to deal with ExpressionSet
#############################################################
removeUnwantedRows<-function(rawData, filters)
{
	cleanRowCount <- 1
	#cleanerData <- rawData #Am assuming this is a copy
	#Get the matrices to work on 
	#featureData
	fdd <-pData(featureData(rawData))
	ed <- assayDataElement(rawData, "exprs")
	fl <- assayDataElement(rawData, "flags")
	cat("Raw Data ", length(fdd[,1]), "\n")
	cat("Number of filters ", length(filters), "\n")
	#Go through for each filter extracting the matrix
	for(x in 1:length(filters))
	{
		#Get the row numbers of all the spots with the filter in the name
		rowsToGo <- grep(filters[x], fdd[, "NAME"], perl=TRUE, ignore.case=FALSE)
		#Get all the rows left 
		allrows <- grep("", fdd[, "NAME"], perl=TRUE, ignore.case=FALSE)
		#Get the set difference i.e. the rows that don't have the matrix in
		#Order matters in setdiff
		rowsToStay <- setdiff(allrows, rowsToGo)
		#Rebuild dataframe without the rowsToGo
		fdd <- fdd[rowsToStay,]
		ed  <- ed[rowsToStay,]
		fl  <- fl[rowsToStay,]
		if(filters[x] == "Empty")
		{
			#Special case for filtering on "" replaced Empty in GAL for 9.2
			cat("Filtering special case empties","\n")
			rowsToGo <- which(fdd[,"NAME"] =="")
			#cat(length(rowsToGo), " empties found\n")
			#Get all the rows left 
			allrows <- grep("", fdd[, "NAME"], perl=TRUE, ignore.case=FALSE)
			#Get the set difference i.e. the rows that don't have the matrix in
			#Order matters in setdiff
			rowsToStay <- setdiff(allrows, rowsToGo)
			#Rebuild dataframe without the rowsToGo
			fdd <- fdd[rowsToStay,]
			ed  <- ed[rowsToStay,]
			fl  <- fl[rowsToStay,]
		}
	}
	#Rebuild the shorter ExpressionSet
	cat("Filtered Data ", length(fdd[,1]), "\n")
	#colnames(fl) = colnames(ed) = pData(phenoData(rawData))
	feats <- new("AnnotatedDataFrame", as.data.frame(fdd))
	varMetadata(feats)<- varMetadata(featureData(rawData))
	cleanerData <- new("ExpressionSet", exprs=ed, featureData=feats, phenoData=phenoData(rawData))
	assayDataElement(cleanerData, "flags") <- fl
	return(cleanerData)
}
#############################################################
#############################################################
#Name:standardRemoveRows
#Description: a standard procedure for removing specific probe rows
#because of the spotting procedure at EMBL, e.g. empty spots, mismatches
#etc.Normally this method is called as a wrapper to removeUnwantedRows
#Author:Jonathon Blake
#Date:30.05.2007
#############################################################
standardRemoveRows <- function(rawData)
{
	filters=c("empty", "Hy3", "_MM1","_MM2", "control", "U6-sn", "No_known","Empty", "Not_designed_for", "No known", "miRPlus", "SNORD", "^-1")
	filteredData <- removeUnwantedRows(rawData, filters)
	return(filteredData)
}
#############################################################
#############################################################
#Name:myForgivingMedian
#Description:Only calculates the median from values that are
#actually present
#Author: Jonathon Blake
#Date:30.11.2007
#############################################################
myForgivingMedian <- function(mat, minSumlength=0)
{
	dmat <- na.omit(mat)
	#A few things have to be done here to prevent it crashing
	#First have to make sure that the length of the matrix is greater than 0 when we take the median
	if(length(dmat)< minSumlength || length(dmat) ==0 )
	{
	#	cat("Array is too short","\n");
		return(NA)
	}
	medVal <- median(dmat)
	#cat("median is ", medVal,  " and sd ", sd(dmat),"\n");
	#Here have to make sure that the length is greater than 1 for an SD
	#I filter here to stop short arrays crashing. It can be that only two readings
	#are present at the start and that you don't lose everything if one is bad.
	#If you are prepared to throw out values that are 1 of two use minSumlength
	if(length(dmat) >1 && abs(medVal) < sd(dmat))
	{
	#	cat("sd is greater than median","\n");
		return(NA)
	}
	#cat("length of medarr = ", length(dmat), "the median is ",median(dmat), "\n")
	return(median(dmat))

}
#############################################################
#############################################################
#Name:naOmitMedian
#Description:Calculates the median from values that are
#actually present Ignores NAs. If madAdjust is true will set
#values to NA is the MAD is greater than the median
#Author: Jonathon Blake
#Date:04.12.2007
#############################################################
naOmitMedian <- function(mat,madAdjust=FALSE)
{
	dmat <- na.omit(mat)
	medmad <- mad(mat, constant=1,na.rm=TRUE) 
	retval <- median(dmat)
	if(is.na(retval))
	{
		return(retval)
	}
	if(madAdjust)
	{
		#cat(retval," ", medmad,"\n")
		if(abs(medmad) > abs(retval))
		{
			#cat("Variation too large median=",retval, "mad=",medmad, "\n")
			retval <- NA;
		}
	}
	#cat(" median=", median(dmat), " mad=",medmad, "\n")
	return(retval)
}
#############################################################
#############################################################
#Name:setIntensityCutoff
#Description: Sets a cutoff value for the intensity. Effectively
#low intensity values have not been trusted. Here anthing below is set
#to NA which effectively says anything below this value is noise.
#Author:Jonathon Blake 
#Date:4.12.2007
#############################################################
setIntensityCutoff <- function(dmat, intensityCutoff)
{
	cat("IntensityCutoff ", intensityCutoff, "\n")
	for(x in 1:length(colnames(dmat)))
	{
		lmat <- dmat[,x] < intensityCutoff
		dmat[lmat,x] <- NA
	}
	return(dmat)
}
#############################################################
#############################################################
#Name:correctForFlags
#Description: Corrects for the flags indicating problems with
#the spot. GENEPIX calls are 0 present and -50, -75, -100 Absent
#This function sets the intensity values associated with these flags
#to NA
#Author:Jonathon Blake
#Date:29.11.2007
#############################################################
correctForFlags<-function(eset, intensityCutoff=0)
{
	#Get the flags matrix
	flag <- assayDataElement(eset, "flags")
	#Get the data that will be updated
	dat <- assayDataElement(eset, "exprs")
	
	#Set the data to NA if flag is less than 0
	for(y in 1:length(colnames(flag)))
	{
		
		badvals <- flag[,y] <0
		dat[badvals, y] <- NA
		if(intensityCutoff >0)
		{
			lmat <- dat[,y] < intensityCutoff
			dat[lmat,y] <- NA
			cat(colnames(dat)[y], " had ", summary(flag[,y]<0)[3], " bad flags  and ", summary(dat[,y]< intensityCutoff)[3], " values below cutoff from ", length(dat[,y])," data points\n")
		}
		else
		{
			cat(colnames(dat)[y], " had ", summary(flag[,y]<0)[3], " bad flags from ", length(dat[,y])," data points\n")
		}
	}
	#if(intensityCutoff >0)
	#{
	#	dat <- setIntensityCutoff(dat, intensityCutoff)
	#}
	correctedData <- new("ExpressionSet", exprs=dat, featureData=featureData(eset), phenoData=phenoData(eset))
	assayDataElement(correctedData, "flags") <- flag
	#return the flag corrected expressionset
	return(correctedData)
}
#############################################################
#############################################################
#Name:summarizeIntensitiesAsMedian
#Description: As the probes are spotted onto the in quaduplet
#or duplicate the values have to be combined in some way. This function
#takes the median of the intensities for the spots. Effectively
#the mean for duplicates. 
#Author:Jonathon Blake
#Date:1.06.2007
#############################################################
summarizeIntensitiesAsMedian <-function(eset,minSumlength=0,madAdjust=FALSE)
{
	#Get the probe IDs
	fd <- pData(featureData(eset))
	#For all the unique IDs for the probes
	uniqueIDs <- unique(fd[,"PROBEID"])
	#Getting the expression data for the calculation
	dat <- assayDataElement(eset, "exprs")

	#For Each ID calculate the median each of the chips
	idrows1 <- grep(uniqueIDs[1], fd[,"PROBEID"])
	idmat1<-dat[idrows1,]
	medrow1 <- apply(idmat1, 2, median)
	cnames <-names(medrow1)
	medmat<- as.matrix(medrow1)
	medmat<-t(medmat)
	tfeats <- idrows1[1]
	
	for(i in 2:length(uniqueIDs))
	{
		#Get the rows for the id and build a matrix to work with
		idrows <- grep(uniqueIDs[i], fd[,"PROBEID"])
		idmat<-dat[idrows,]
		#medrow <- apply(idmat, 2, naOmitMedian)
		#medrow <- apply(idmat, 2, myForgivingMedian, minSumlength=minSumlength)
				#
		medrow <- apply(idmat, 2, naOmitMedian, madAdjust)
		#cat(uniqueIDs[i], " ") 
		#cat("Binding row together", "\n")
		medmat<-rbind(medmat,medrow)	
		tfeats <- c(tfeats, idrows[1])
	}
	#Rebuild the names and colnames on the dataframe
	colnames(medmat) <- cnames
	row.names(medmat) <-uniqueIDs
	tfs <- fd[tfeats,2:3]
	colnames(tfs) <- cbind("PROBEID", "NAME")
	row.names(tfs) <-uniqueIDs
	tfd <- new("AnnotatedDataFrame", as.data.frame(tfs))
	#labdata <- as.data.frame(c("Exiqon Probe ID", "mirTarget Name"))
	labdata <- as.data.frame(c("Exiqon Probe ID", "miRNA Name"))
	colnames(labdata) <- 'labelDescription'
	varMetadata(tfd)<- labdata
	summedSet <- new("ExpressionSet", exprs=medmat, phenoData= phenoData(eset), featureData=tfd)
	return(summedSet)
	
}
#############################################################
#NORMALIZATION FUNCTIONS
#############################################################
#############################################################
#Name:normalizePerChipMedian
#Description:Normalizes each chip to the median of the chip
#Takes an expression set and normalizes each column to its median
#Author: Jonathon Blake
#Date: 30.11.2007
#############################################################
normalizePerChipMedian <-function(eset)
{
	dmat <- assayDataElement(eset, "exprs")
	#chipMeds <- apply(dmat, 2, naOmitMedian)
	medmat<- NULL
	#for(x in 1:length(chipMeds))
	for(x in 1:length(colnames(dmat)))
	{
		medt <- na.omit(dmat[,x])
		mednormed <- dmat[,x]/median(medt)
		#cat("The median is ", median(medt)," the length is  ", length(medt),"\n")
		#mednormed <- dmat[,x]/chipMeds[x]
		
		medmat<- cbind(medmat, mednormed)
	}
	colnames(medmat) <- colnames(dmat)	
	row.names(medmat) <-row.names(dmat)
	normedSet <- new("ExpressionSet", exprs=medmat, phenoData= phenoData(eset), featureData=featureData(eset))
	return(normedSet)
}

#############################################################
#PLOTTING FUNCTIONS
#############################################################
#Name: boxPlotData
#Description: Takes a data matrix normally the expression matrix
#from an ExpressionSet and produces a boxplot. It omits nas which have
#caused a number of problems in the boxplots in the past.
#Author: Jonathon Blake
#Date: 30.11.2007
#0.5.1 used the rainbow function for generating the colour arrays
#############################################################
boxplotData <-function(dmat, exptname, dlevel)
{

	filedes <- paste(exptname, dlevel, sep="_")
	jpgname <- paste(filedes, ".jpg", sep="")
	jpeg(filename=jpgname, height=600, width=800)
#0.5.1

boxplot(as.data.frame(log2(dmat)), xlab="Samples", ylab="Log2 Intensity", main=filedes, col=rainbow(length(colnames(dmat))))
	dev.off()
	cat("File saved to ", file.path(getwd(), jpgname),"\n", sep="")
	gc()
	
}
#############################################################
#############################################################
#Name: boxPlotDataNoFile
#Description: Takes a data matrix normally the expression matrix
#from an ExpressionSet and produces a boxplot. It omits nas which have
#caused a number of problems in the boxplots in the past.
#Author: Jonathon Blake
#Date: 30.11.2007
#0.5.1 used the rainbow function for generating the colour arrays
#############################################################
boxplotDataNoFile <-function(dmat, exptname, dlevel)
{
	filedes <- paste(exptname, dlevel, sep="_")
	boxplot(as.data.frame(log2(dmat)), xlab="Samples", ylab="Log2 Intensity", main=filedes, col=rainbow(length(colnames(dmat))))
	gc()
	
}
#############################################################
#Name: panelCor
#Description: Adds a pearson correlation value to the scatter plots
#The code is taken from the panel manpage and modified
#Author:Jonathon Blake
#Date:4.02.08
#############################################################
panelCor <-function(x,y,digits=2, prefix="r=")
{
	ptest <-cor.test(x, y, na.action=na.omit)
	ptest
 	r <- ptest$estimate;
	txt <- format(c(r, 0.123456789), digits=digits)[1];
	txt <- paste(prefix, txt, sep="");
    	cex <- 0.8/strwidth(txt)
	#text(0.5, 0.5, txt, cex = cex * r)
	#text(naOmitMedian(x), naOmitMedian(y), txt, cex)
	text(2, 2, txt, cex)
}
#############################################################
#Name: plotIntensitiesScatter
#Description: Works to plot intensities of probes pairwise
#against each other 
#Author:Jonathon Blake
#Date: 11.07.07
#############################################################
plotIntensitiesScatter <-function(dmat, controls=NULL, exptname, maintitle)
{
	iname <- paste(exptname, maintitle, sep="\n")
	matlen <- length(colnames(dmat))
	
	#cat( "length of names(dmat) is ", matlen)
	#cat( "length of dmat is ", length(dmat))
	jname <-paste(exptname, maintitle, sep="_")
	jpgname <-paste(jname, ".jpg", sep="")
	jpeg(filename=jpgname, height=1000, width=1000)
	pairs(log10(dmat), panel=function(...) {points(..., pch=".", asp=1); abline(a=0, b=1, col="red");  panelCor(...)})
	dev.off()
	cat("File saved to ", file.path(getwd(), jpgname),"\n", sep="")


}
#############################################################
#############################################################
#Name: workedExampleMedianNormalize
#Description: A demonstration of how the code can be used
#Author: Jonathon Blake
#Date: 30.05.2007
#############################################################
workedExampleMedianNormalize <-function(exptname, intensityCutoff=0, datadir=".",minSumlength=0, madAdjust=FALSE)
{
	#Get the data out of the genepix files and into a dataframe
	#The call below is searching the current directory for genepix files 
	myRawData <- parseRawData(datadir)
	
	#################
	#Filter the data frame to remove whatever rows have controls, are empty etc
	#################
	filteredData <-standardRemoveRows(myRawData)
	outputAnnotatedDataMatrix(filteredData, exptname, "FilteredData", "exprs")
	
	#################
	#Correct for flags and remove the flag columns
	#################
	flagcorData <- correctForFlags(filteredData, intensityCutoff)
	boxplotData(exprs(flagcorData), exptname, "FlagCorrected")
	outputAnnotatedDataMatrix(flagcorData, exptname, "FlagCorrectedIntensity", "exprs")
	
	#################
	#Summarize the data for each of the 4 spots on the chip into medians
	#################
	summedData <-summarizeIntensitiesAsMedian(flagcorData, minSumlength,madAdjust)
	boxplotData(exprs(summedData), exptname, "Summarized")
	outputAnnotatedDataMatrix(summedData, exptname, "summarizedIntensity", "exprs")
	plotIntensitiesScatter(exprs(summedData),NULL, exptname,"SummarizedScatter")
	
	#################
	#Get the median normalized data
	#################
	mednormedData <- normalizePerChipMedian(summedData)
	boxplotData(exprs(mednormedData), exptname, "Median Normalized")
	plotIntensitiesScatter(exprs(mednormedData),NULL, exptname,"MedNormedScatter")
	outputAnnotatedDataMatrix(mednormedData, exptname, "medNormedIntensity", "exprs")
	
	#################
	#Return normalized data
	#################
	return(mednormedData)

}
#############################################################
#############################################################
#Name: workedExampleNotNormalizedData
#Description: A demonstration of how the code can be used
#Author: Jonathon Blake
#Date: 30.05.2007
#############################################################
workedExampleNotNormalizedData <-function(exptname, intensityCutoff=0, datadir=".",minSumlength=0, madAdjust=FALSE)
{
	#Get the data out of the genepix files and into a dataframe
	#The call below is searching the current directory for genepix files 
	myRawData <- parseRawData(datadir)
	
	#################
	#Filter the data frame to remove whatever rows have controls, are empty etc
	#################
	filteredData <-standardRemoveRows(myRawData)
	outputAnnotatedDataMatrix(filteredData, exptname, "FilteredData", "exprs")
	
	#################
	#Correct for flags and remove the flag columns
	#################
	flagcorData <- correctForFlags(filteredData, intensityCutoff)
	boxplotData(exprs(flagcorData), exptname, "FlagCorrected")
	outputAnnotatedDataMatrix(flagcorData, exptname, "FlagCorrectedIntensity", "exprs")
	
	#################
	#Summarize the data for each of the 4 spots on the chip into medians
	#################
	summedData <-summarizeIntensitiesAsMedian(flagcorData,minSumlength, madAdjust)
	boxplotData(exprs(summedData), exptname, "Summarized")
	outputAnnotatedDataMatrix(summedData, exptname, "summarizedIntensity", "exprs")
	plotIntensitiesScatter(exprs(summedData),NULL, exptname,"SummarizedScatter")
	
	
	#################
	#Return summarized data for further normalization
	#################
	return(summedData)

}
#############################################################


