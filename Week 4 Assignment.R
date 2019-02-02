#Run this script by the individual steps, one at a time; else it may crash.

#install and load packages
install.packages("BiocManager")
library(BiocManager)
install(c("sangerseqR", "annotate", "genbankr"))
library(sangerseqR)
install.packages("seqinr")
library(seqinr)

# 1. Extract and call the squence from the .ab1 file provided

ITS = read.abif("./Data/DNA_Barcoding/1Ipos_F_P1815443_064.ab1") # Read
ITSseq <- sangerseq(ITS) # Extract
SeqX<-makeBaseCalls(ITSseq) # Call 

# 2. Use regular expressions to slice, paste and extract the ‘Primary sequence’ only.

PrimSeq1=SeqX@primarySeq #slice and extract
PrimSeq1 #paste

# 3. Convert the sequence to a FASTA format, in which the first line starts with > and contains the file name and the second line contains the primary sequence.

write.fasta(sequences = PrimSeq1, names = "1Ipos_F_P1815443_064.ab1", nbchar = 80, file.out = "1Ipos_F_P1815443_064.fasta")
cat(paste("> 1Ipos_F_P1815443_064.ab1","\n", read.fasta(file = "1Ipos_F_P1815443_064.fasta", as.string = TRUE, forceDNAtolower = FALSE, seqonly = T, strip.desc = FALSE, apply.mask = TRUE)))

# 4. Loop through every file in the Data folder, do steps 1-3, and save the output as a vector of strings in the FASTA format.
allab1s=list.files("./Data/DNA_Barcoding")
for(i in allab1s){
  #ITS = read.abif(paste("./Data/DNA_Barcoding/",i, sep = ""))
  #ITSseq <- sangerseq(ITS) # Extract
  ITSseq = readsangerseq(paste("./Data/DNA_Barcoding/",i, sep = ""))
  SeqX<-makeBaseCalls(ITSseq) # Call
  PrimSeq1=SeqX@primarySeq #slice and extract
  PrimSeq1 #paste
  fastafile="fastafile"
  write.fasta(sequences = PrimSeq1, names = i, nbchar = 80, file.out = fastafile, open="a")
  cat(paste(">", i,"\n", sep="", read.fasta(file = fastafile, as.string = TRUE, forceDNAtolower = FALSE, seqonly = T, strip.desc = FALSE, apply.mask = TRUE),"\n\n"))
}

# 5. Save the BarcodePlateStats.csv to your Data folder and use it to exclude from your loop any sequences that do not pass the quality check.

library(dplyr) #install dplyr
qual=read.csv("./Data/BarcodePlateStats.csv") #read in BarcodePlateStats .csv file
names(qual)
qualtrue=subset(qual, Ok == "TRUE") #create subset of data with Ok = TRUE
for(i in allab1s){
  if(i %in% qualtrue$Chromatogram){ #exclude files with FALSE condition under "Ok"
    ITSseq = readsangerseq(paste("./Data/DNA_Barcoding/",i, sep = ""))
    SeqX<-makeBaseCalls(ITSseq) # Call
    PrimSeq1=SeqX@primarySeq #slice and extract
    PrimSeq1 #paste
    fastafile="fastafile"
    write.fasta(sequences = PrimSeq1, names = i, nbchar = 80, file.out = fastafile, open="a")
    cat(paste(">", i,"\n", sep="", read.fasta(file = fastafile, as.string = TRUE, forceDNAtolower = FALSE, seqonly = T, strip.desc = FALSE, apply.mask = TRUE),"\n\n"))
  }
}
