#=============================================================================================================#
# Script created by Mark Christie, all rights reserved, contact at markchristie1500@gmail.com
# Script created in version R 4.1.0 on 07/292021
# Script modified by Andy Lee in version R ______ on Oct. 16, 2022. for kellet's whelk data set
# This script: merges blastx output from swissprot and nr databases
# Usage notes: 
#============================================================================================================#
# Set working directory, import packages, source functions, initialize global variables
library('here')
here::i_am("merge.R")
here()

list.files()

sprot <- read.delim("headswissprot.txt", header = FALSE) # read delim is important for this irregular file format
nrow(sprot)

nr    <- read.delim("headnr.txt", header = FALSE) # read delim is important for this irregular file format
nrow(nr)

#write.table(offs, "temp".txt", col.names = TRUE, sep="\t", append = FALSE)
#=============================================================================================================#

# add column indicating database source
head(sprot)
sprot <- cbind(sprot, "swissprot")

# add column indicating database source
head(nr)
nr <- cbind(nr, "nonredundant")

# add column names
cnames          <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "ssciname", "scomname", "stitle", "database")
colnames(sprot) <- cnames
colnames(nr)    <- cnames

dat <- rbind(sprot, nr)

## create unique unmerical index for contig names 
head(dat)
col1 <- dat[, 1]
col2 <- as.numeric(factor(col1))

contigs <- 1:length(unique(col2)) # there are 779 DEGs and single FASTA files
cols    <- unique(col2)

## # calculate number and identifies of missing fastas (i.e., fastas for which there was no hit in swiss prot or nr)
## m1      <- match(contigs, cols)
## missing <- which(is.na(m1) == TRUE)
## missing
## length(missing)
## write.table(missing, "unannotated_contigs.txt", col.names = TRUE, sep="\t", append = FALSE)

# add fasta number to all rows
dat2 <- cbind(col2, dat)
table(dat2[, 1])

# write table for all output
#write.table(dat2, "all_annotations.txt", col.names = TRUE, sep="\t", append = FALSE) # all annotations (still filtered by e-value set during blastx run)

# write table for "best hits"
# "best hits" are as follows:
# 1. lowest evalue
# 2. highest bit score (useful for standardizing across databases, here swiss prot and ncbi)
# 3. (genes$pident * .01)*genes$length   == percent identity (percent of target sequence that matches) * length of target sequence =  total length of target sequence that matches query
# head(dat2)
# contigs <- unique(dat2[, 1])
# OUT <- NULL
# for(n in contigs){
#   genes <- dat2[dat2[, 1] == n, ]
#   gen1  <- genes[which.min(genes$evalue), ]
#   gen2  <- genes[which.max(genes$bitscore), ]
#   ptarg <- (genes$pident * .01)*genes$length 
#   gen3  <- genes[which.max(ptarg), ]
#   output <- rbind(gen1, gen2, gen3)
#   OUT <- rbind(OUT, output)
#   }

#write.table(OUT, "top_annotations_final.txt", col.names = TRUE, sep="\t", append = FALSE) # top annotations 

# select top e-value only
head(dat2)
contigs <- unique(dat2[, 1])
OUT <- NULL
for(n in contigs){
  genes <- dat2[dat2[, 1] == n, ]
  gen1  <- genes[which.min(genes$evalue), ]
  output <- gen1
  OUT <- rbind(OUT, output)
}

#sort
dat <- OUT
dat <- dat[order(dat[, 1]), ]

# Add in all gene IDS
contigs <- 1:length(unique(col2))

OUT2 <- NULL
for(g in contigs){
  id   <- g
  m1   <- match(g, as.numeric(factor(dat[, 1])))
  out  <- dat[m1, ]
  OUT2 <- rbind(OUT2, out)
}

nrow(OUT2)
OUT2 <- cbind(contigs, OUT2)
write.table(OUT2, "top_e_value_hit.txt", col.names = TRUE, sep="\t", append = FALSE) # top annotations 

#================#
