# install.packages("optparse")
# library(optparse)
library(stringr)
library(reshape)

#########################################################################################
# File Constants
#########################################################################################

fileOrths <- "/auto/sahara/namib/home/mprice/data/mutagenesis/MR1/MR1_Shewanella_orths"
fileMR1Genes <- "/auto/sahara/namib/home/mprice/data/FEBA/g/MR1/genes.tab"
fileMR1 <- "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//MR1_Results.tab"
fileEssentialEColi <- "/usr2/people/mprice/tmp/shg"

l_strFiles <- list("/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//SB2B_Results.tab",
                   "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//ANA3_Results.tab",
                   "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//PV4_Results.tab")


#########################################################################################
# Functions
#########################################################################################

getStrain<- function(strName){
  # Gets the strain name from its results file.
  strStrain <- str_match(strName, ".*insertion_counts//(.*)_Results.tab")[2]
  return(strStrain)
}

prepFileForMerge<- function(fileCounts,dfOrths){
  # Converts the locus on file to correspond to MR1.
  
  strStrain <- getStrain(fileCounts)
  # Cut down the main file to the two appropriate columns
  dfOrths_Small <- subset(dfOrths, select=c("MR1",strStrain))
  
  # Read in the insertions
  dfInserts <- read.table(fileCounts,header=FALSE)
  colnames(dfInserts)<- c(strStrain,"Insertions")
  
  # Merge on the MR1 locus
  dfMerged <- merge(dfInserts,dfOrths_Small,by=strStrain)
  dfMerged <- subset(dfMerged,select=c("MR1","Insertions"))
  dfMerged$Insertions <- as.numeric(as.character(dfMerged$Insertions))
  colnames(dfMerged) <- c("MR1_locus",strStrain)
  
  
  return(dfMerged)
}

IsZeroOrNA<- function(x){
  if(x==0 | is.na(x)){
    return(TRUE)}
  else{
    return(FALSE)}
}

countZeroNA<- function(x){
  sum(sapply(x,IsZeroOrNA))
}

#######################################################


###########################################################################################
# Load and format file of orthologues.

# "c_vstrMR1Header" obtained from 1/9/2014 email from Morgan. It corresponds to the file at
# ~mprice/data/mutagenesis/MR1/MR1_Shewanella_orths.
# If we need the strain names for the other columnbs, they can be obtained from 
# http://www.ncbi.nlm.nih.gov/taxonomy

c_vstrMR1Header <- c("MR1","orth60480.locusId","orth60481.locusId","orth94122.locusId","orth225849.locusId","orth318161.locusId","orth318167.locusId","orth319224.locusId","orth323850.locusId","orth325240.locusId","orth326297.locusId","orth351745.locusId","orth392500.locusId","orth398579.locusId","orth399599.locusId","orth402882.locusId","orth407976.locusId","orth425104.locusId","orth458817.locusId","orth637905.locusId")
vstrOrthHeader <- str_replace(c_vstrMR1Header, "orth(.*)\\.locusId", "\\1")
vstrOrthHeader <- str_replace(vstrOrthHeader, "326297", "SB2B")
vstrOrthHeader <- str_replace(vstrOrthHeader, "94122", "ANA3")
vstrOrthHeader <- str_replace(vstrOrthHeader, "323850", "PV4")



dfOrths <- read.table(fileOrths, header=FALSE, sep = "\t",quote="")
colnames(dfOrths) <- vstrOrthHeader

###########################################################################################
# Load in file of MR1 Gene information.

dfMR1Genes <- read.table(fileMR1Genes,sep="\t",header=TRUE,quote="")
colnames(dfMR1Genes)[1] <- "MR1_locus"
write.table(dfMR1Genes, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/MR1_Check.tab",sep="\t",row.names=FALSE)
print(paste0("There are ",nrow(dfMR1Genes)," MR1 genes in the main file."))


###########################################################################################
# Merge the files.

dfMR1 <- read.table(fileMR1,sep="\t",header=TRUE)
colnames(dfMR1)<- c("MR1_locus","MR1")

l_dfStrains <- lapply(l_strFiles, prepFileForMerge,dfOrths=dfOrths)

dfOrths_Small <- subset(dfOrths, select=c("MR1"))
colnames(dfOrths_Small) <- c("MR1_locus")

dfMR1Orths <- merge(dfOrths_Small,dfMR1,all.x=TRUE,by="MR1_locus")

l_dfToMerge <- c(list(dfMR1Orths),l_dfStrains)

# Code modified from:
# http://stackoverflow.com/questions/2209258/merge-several-data-frames-into-one-data-frame-with-a-loop
dfAll <- Reduce(function(x, y) merge(x, y, all.x=T, 
          by="MR1_locus"), l_dfToMerge, accumulate=F)


dfAll <- merge(dfAll,dfMR1Genes,all.x=TRUE,by="MR1_locus")
dfAll <- dfAll[with(dfAll, order(MR1,SB2B,ANA3,PV4)), ]

write.table(dfAll, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/All.tab",sep="\t",row.names=FALSE)

############################################################################################
dfEColi <- read.table(fileEssentialEColi,sep="\t",header=TRUE,quote="")
print(paste0("There are ",nrow(dfEColi), " rows in the EColi dataframe."))

dfStrainsEColi <- merge(dfAll,dfEColi,by.x="sysName",by.y="sysName2",all.x=TRUE)

dfStrainsEColi$MorganEssential <- (dfStrainsEColi$Essential=="E" | dfStrainsEColi$Class==1)
dfStrainsEColi <- dfStrainsEColi[dfStrainsEColi$MorganEssential==FALSE & !is.na(dfStrainsEColi$MorganEssential),]

dfStrainsEColi$WithoutInserts <- apply((dfStrainsEColi[,3:6]),1,countZeroNA)


print(dfStrainsEColi[2,3:6])
print(countZeroNA(dfStrainsEColi[2,3:6]))
print(countZeroNA(c(0,0,0,0)))




write.table(dfStrainsEColi,"/auto/sahara/namib/home/jkaminski/feba_essential/tmp/EColi.tab",sep="\t",row.names=FALSE )


# These are genes that are not essential in ecoli, but appear to be essential in our genomes.
dfGenesOfInterest <- dfStrainsEColi[dfStrainsEColi$WithoutInserts==4 ,]


write.table(dfGenesOfInterest,"/auto/sahara/namib/home/jkaminski/feba_essential/tmp/EssentialShewanella_Not_EColi.tab",sep="\t",row.names=FALSE )
