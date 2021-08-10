# Title: Data Mining and Subsampling in Parallel Evolution in Influenza Viral RNA Sequences
# AUTHORS

# Jennifer E. Jones <jennyjones_xyz@pitt.edu>
# Erik S. Wright <eswright@pitt.edu>
# Seema S. Lakdawala <lakdawala@pitt.edu>
  
# This analysis requires the DECIPER package. This package and its dependencies
# can be downloaded from Bioconductor at http://bioconductor.org/packages/release/bioc/html/DECIPHER.html.

# Load libraries.

library(DECIPHER)

sessionInfo() # Need DECIPHER >= v2.20.0

# Set working directory.

setwd("<<PATH TO Parallel-Evolution>>/") # Needs trailing slash.
# Example "~/Desktop/Virus/flu/Parallel-Evolution-main/"

# FASTA FILE IMPORT AND CLEANUP

# Load FASTA sequences for analysis. Datasets analyzed in this manuscript can be found at 
# https://github.com/Lakdawala-Lab/Parallel-Evolution-Of-Influenza-Viral-RNA/Preprocessed FASTA Files/.
# Analysis of the H3N2 virus sequences are presented as an example here (available in the sub-folder 'H3N2 Virus Sequences'.)

subtype <- "H3N2"
segments <- c("1-PB2", "2-PB1", "3-PA", "4-HA", 
"5-NP", "6-NA", "7-M", "8-NS")
nSegments <- length(segments)

files <- paste("./Human influenza", subtype, "segment", segments, "sequences.fasta.gz")

# FASTA file QC.

length <- length(files)
mRNA <- vector(mode = "list", length = length)
vRNA <- vector(mode = "list", length = length)

for (i in seq_along(vRNA)) {
      mRNA[[i]] <- readDNAStringSet(files[i]) # IAV fasta files are downloaded as positive-sense (mRNA) sequences. 
      vRNA[[i]] <- reverseComplement(mRNA[[i]]) # Convert sequence files to the genomic (vRNA) sequences.
      vRNA[[i]] <- vRNA[[i]][grep(paste0("|Subtype:", subtype, "|"), names(vRNA[[i]]), fixed=TRUE)] # Ensure that all fasta files contain the desired subtype (H1N1 or H3N2).
      vRNA[[i]] <- vRNA[[i]][grep("Host:Human", names(vRNA[[i]]), fixed=TRUE)] # Ensure that all fasta files contain human origin IAV.
      
      # Extract the strain names from the FASTA files. This step is required to ensure accurate concatenation of full-length genomes.
      names(vRNA[[i]]) <- gsub(".+\\|Strain Name:(.+?)\\|Segment.+", "\\1", names(vRNA[[i]])) 
      
      # Remove duplicate entries from FASTA files.
      d <- which(!duplicated(names(vRNA[[i]])))
      vRNA[[i]] <- vRNA[[i]][d]
}

# Remove FASTA files with incomplete genomic sequences from analysis.

t <- table(unlist(sapply(vRNA, function(x) unique(names(x)))))
t <- names(t)[t == nSegments]
vRNA <- lapply(vRNA, `[`, t)

# Assemble full-length genomes by concatenation.

vRNA[[9]] <- do.call(xscat, vRNA)
names(vRNA[[9]]) <- names(vRNA[[1]])
out_files <- c(gsub("sequences", "vRNA sequences", files),
	gsub("segment .+ sequences", "concatenated full-length vRNA sequences", files[1]))

# Optionally, write all files to the designated working directory as a potential stopping point.

for (i in seq_along(vRNA))
	writeXStringSet(vRNA[[i]], file =  out_files[i], compress=TRUE)
# Reading files back in after stopping. 
length <- length(out_files)
vRNA <- vector(mode = "list", length = length)
vRNA <- lapply(out_files, readDNAStringSet)

# Selecting strains from specific time periods for further analysis.

year <- as.numeric(gsub(".+/([0-9])", "\\1", names(vRNA[[1]])))
w <- which(year > 1994 & year < 2005) # 1995 - 2004
vRNA1 <- lapply(vRNA, `[`, w)
w <- which(year > 2004 & year < 2015) # 2005 - 2014
vRNA2 <- lapply(vRNA, `[`, w)

vRNA1 <- lapply(vRNA1, AlignSeqs, processors = NULL)
vRNA2 <- lapply(vRNA2, AlignSeqs, processors = NULL)

# CLUSTERING INTO OPERATIONAL TAXONOMIC UNITS (OTUs) AND SEQUENCE SELECTION

# Calculate distances between full-length concatenated sequences 
# and compare clusters with different cutoffs ranging from 95-99% sequence identity.

d1 <- DistanceMatrix(vRNA1[[9]], type="dist", correction="JC", processors = NULL)
d2 <- DistanceMatrix(vRNA2[[9]], type="dist", correction="JC", processors = NULL)

otu1 <- IdClusters(d1, method = "NJ", cutoff = c(0.01, 0.02, 0.03, 0.04, 0.05), type = "clusters", myXStringSet = vRNA1[[9]], processors = NULL)
otu2 <- IdClusters(d2, method = "NJ", cutoff = c(0.01, 0.02, 0.03, 0.04, 0.05), type = "clusters", myXStringSet = vRNA2[[9]], processors = NULL)

# View the number of clusters in each species tree with each cutoff.

sapply(otu1, max)
sapply(otu2, max)

# Choose the desired cutoff (in this case, 97% sequence identity was selected).

otu1 <- otu1[, 3, drop=FALSE]
otu2 <- otu2[, 3, drop=FALSE]

# Write clustering data to working directory.

write.table(otu1, file = "H3N2 1995-2004 Clusters.csv", sep=",")
write.table(otu2, file = "H3N2 2005-2014 Clusters.csv", sep=",")

# Manually inspect sequence quality and choose representative sequences from each cluster. Then, subset the alignments 
# of the sequences selected. The sequences analyzed here are indicated in the 'H3N2 1995-2004 Strains Analyzed.csv' and 
# 'H3N2 2005-2014 Strains Analyzed.csv', which can be found in the following sub-folder:
# https://github.com/Lakdawala-Lab/Parallel-Evolution-Of-Influenza-Viral-RNA/Preprocessed FASTA Files/H3N2 Virus Sequences/.

id <- read.csv("H3N2 1995-2004 Strains Analyzed.csv")
id <- id[2:length(id)]

for (i in seq_along(files)) {
  for (j in seq_along(id)) {
    m <- match(names(vRNA1[[i]]), id[, j])
    w <- which(!is.na(m))
    writeXStringSet(vRNA1[[i]][w], file= paste0(substring(files[i], 1, nchar(files[i]) - nchar("sequences.fasta.gz")), "vRNA MSA - 1995-2004 OTU ", j, ".fasta.gz", sep=""))
  }
}

id <- read.csv("H3N2 2005-2014 Strains Analyzed.csv")
id <- id[2:length(id)]

for (i in seq_along(files)) {
  for (j in seq_along(id)) {
    m <- match(names(vRNA2[[i]]), id[, j])
    w <- which(!is.na(m))
    writeXStringSet(vRNA2[[i]][w], file= paste0(substring(files[i], 1, nchar(files[i]) - nchar(".fasta.gz")), " OTU ", j, ".fasta.gz", sep=""))
  }
}

# Fully processed FASTA files are now ready for tree reconstruction and analysis of tree similarity. Source code to complete this analysis
# is provided in the 'Analysis of Tree Similarity' sub-folder.





