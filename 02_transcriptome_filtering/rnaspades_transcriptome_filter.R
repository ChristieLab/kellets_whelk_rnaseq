#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script modified by Andy Lee to be run on beagle and test different filtering parameters 
# Script created in version R 4.1.0 on xx/xx/2021
# This script: Takes a fasta file (e.g., a transcriptome) and reports average sequence length and allows for removing sequences < n bp
# Usage notes: Must set a cut-off value 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
library('seqinr')
library('plyr')
library('dplyr')

fasta <- read.fasta("rnaspades_assembly_calpoly.fasta", as.string = TRUE)
head(fasta)

OUT <- NULL
for(n in 1:length(fasta)){
  fasta1       <- fasta[[n]][1]
  fasta.length <- nchar(fasta1)  
  bases        <- unlist(strsplit(fasta1, ""))
  gcat         <- table(bases)                           
  out          <- as.data.frame(cbind(n, fasta.length, t(gcat)))  
  # if(length(gcat) < 4) {gcat <- c(0,0,0,0)} #check to see how many contigs are missing gcat values
  OUT          <- rbind.fill(OUT, out) #rbind.fill adds NA to missing columns 
}

# gc content
gc <- (OUT[, 4] + OUT[, 5])/(OUT[, 3] + OUT[, 4] + OUT[, 5] + OUT[, 6])

# save some plots
pdf("rnaspades_calpoly.pdf", height=8, width=16) # open a blank pdf with specification 
hist(OUT[, 2], breaks = 100, xlim =c(0,1000), xlab = "contig size") 
hist(gc, breaks = 50,xlim =c(0,1000), xlab ="gc content")   # consider filtering out on this criteria too
plot(OUT[, 2], gc,xlim =c(0,1000))  # correlation between gc content and contic length???
dev.off() # close the pdf 

OUT <- cbind(OUT,gc) # add gc content to OUT dataframe 

# create reduced fasta file
#filter 

#keep_l400    <- which(OUT[, 2] >= 400) # filter out length < 400
#keep_gc30   <- which(OUT[, 7] >= .30) # 30 percent gc content
#keep_l1000    <- which(OUT[, 2] >= 1000) # filter out length < 1000

#fasta_l400 <- fasta[keep_l400]
#fasta_gc30 <- fasta[keep_gc30]
#fasta_l1000 <-fasta[keep_l1000]

#write.fasta(sequences=fasta_l400,names=names(fasta_l400),file.out= "rnaspades_assembly_calpoly_filtered_l400.fasta", as.string = FALSE) # may need to set as.string to true; also note previous capitalized letters now lowercase
#write.fasta(sequences=fasta_gc30,names=names(fasta_gc30),file.out= "rnaspades_assembly_calpoly_filtered_gc30.fasta", as.string = FALSE)
#write.fasta(sequences=fasta_l1000,names=names(fasta_l1000),file.out= "rnaspades_assembly_calpoly_filtered_l1000.fasta", as.string = FALSE)
write.csv(OUT, file="rnaspades_assembly_dataframe.csv")
