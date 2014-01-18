library(stringr)

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
  colnames(dfMerged) <- c("MR1_locus",strStrain)
  
  
  return(dfMerged)
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

fileOrths <- "/auto/sahara/namib/home/mprice/data/mutagenesis/MR1/MR1_Shewanella_orths"

dfOrths <- read.table(fileOrths, header=FALSE, sep = "\t",quote="")
colnames(dfOrths) <- vstrOrthHeader

###########################################################################################
# Load in file of MR1 Gene information.

dfMR1Genes <- read.table("/auto/sahara/namib/home/mprice/data/FEBA/g/MR1/genes.tab",sep="\t",header=TRUE,quote="")
colnames(dfMR1Genes)[1] <- "MR1_locus"
write.table(dfMR1Genes, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/MR1_Check.tab",sep="\t",row.names=FALSE)
print(paste0("There are ",nrow(dfMR1Genes)," MR1 genes in the main file."))


###########################################################################################
# Merge the files.
# (Change this part to use "lapply".
dfSB2B <- prepFileForMerge("/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//SB2B_Results.tab",dfOrths)
dfANA3 <- prepFileForMerge("/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//ANA3_Results.tab",dfOrths)
dfMR1 <- read.table("/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts//MR1_Results.tab",sep="\t",header=TRUE)
colnames(dfMR1)<- c("MR1_locus","MR1")



dfOrths_Small <- subset(dfOrths, select=c("MR1"))
colnames(dfOrths_Small) <- c("MR1_locus")



dfAll <- merge(dfOrths_Small,dfMR1,all.x=TRUE,by="MR1_locus")
dfAll <- merge(dfAll,dfSB2B,all.x=TRUE,by="MR1_locus")
dfAll <- merge(dfAll,dfANA3,all.x=TRUE,by="MR1_locus")
dfAll <- merge(dfAll,dfMR1Genes,all.x=TRUE,by="MR1_locus")
print(paste0("There are ",nrow(dfAll)," genes in the final results file."))

dfAll <- dfAll[with(dfAll, order(MR1,SB2B,ANA3)), ]





write.table(dfSB2B, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/SB2B.tab",sep="\t",row.names=FALSE)
write.table(dfANA3, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/ANA3.tab",sep="\t",row.names=FALSE)
write.table(dfMR1, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/MR1.tab",sep="\t",row.names=FALSE)
write.table(dfAll, file = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/All.tab",sep="\t",row.names=FALSE)



# # Load in Insertion Counts
# dirInserts = "/auto/sahara/namib/home/jkaminski/feba_essential/tmp/insertion_counts/"
# 
# lFiles=list.files(path=dirInserts, full.names=TRUE)
# print(lFiles)
# 
# lstrStrains <- lapply(lFiles,getStrain)
# 
# ldfInserts <- lapply(lFiles, read.table,header=TRUE)
# 
# print(str(ldfInserts))
# print(lstrStrains)
# 
# # Merge MR1 id onto each df in ldfInserts.
# # Merge all of these together by the MR1 id.
# 
# #ldfInsertsMR1 <- lapply(ldfInserts, merge,header=TRUE)


