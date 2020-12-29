library(maanova)

ghanaDDT <- read.madata(datafile="ghanaDDT.data.txt", 
	designfile="DDTdesignfile.txt", 
	arrayType="twoColor", header=TRUE, probeid=1, genename=2, 
	systematicname=3, description=4, status=5, sequence=6, metarow=7, 
	metacol=8, row=9, col=10, intensity=11, spotflag=TRUE, 
	avgreps=0, log.trans=TRUE)

setwd ("C:/Martin/Data/Craig/Tororo microarray")

tororo <- read.madata(datafile="Tororo.data.txt", 
	designfile="Tororodesignfile.txt", arrayType="twoColor", 
	header=TRUE, probeid=1, genename=2, systematicname=3, 
	description=4, status=5, sequence=6, metarow=7, metacol=8, 
	row=9, col=10, intensity=11, spotflag=TRUE, avgreps=0, log.trans=TRUE)

riplot(tororo)
