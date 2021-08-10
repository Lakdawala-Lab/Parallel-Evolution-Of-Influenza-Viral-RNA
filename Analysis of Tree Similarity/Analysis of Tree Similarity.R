# Title: Analysis of Tree Similarity in Parallel Evolution in Influenza Viral RNA SequencesPre
# AUTHORS

# Jennifer E. Jones <jennyjones_xyz@pitt.edu>
# Erik S. Wright <eswright@pitt.edu>
# Seema S. Lakdawala <lakdawala@pitt.edu>

# This analysis requires the DECIPER package. This package and its dependencies
# can be downloaded from Bioconductor at http://bioconductor.org/packages/release/bioc/html/DECIPHER.html.

# Load libraries.

library(DECIPHER)
library(ape)
library(phangorn)
library(RColorBrewer)
library(lattice)
library(plotrix)
library(TreeDist)

sessionInfo() # Need DECIPHER >= v2.18.1

# Set working directory.

setwd("<<PATH TO Parallel-Evolution>>/") # Needs trailing slash.
# Example "~/Desktop/Virus/flu/Parallel-Evolution-main/"

# RECONSTRUCTING PHYLOGENETIC TREES

# Analysis of one set of H3N2 virus trees is illustrated here.

# First perform model testing on processed FASTA files for maximum-likelihood method. FASTA files can be loaded from
# https://github.com/Lakdawala-Lab/Parallel-Evolution-Of-Influenza-Viral-RNA/Processed Fasta Files/H3N2 Virus Sequences/

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 1.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 2.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 2.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 3.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 3.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 4.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 4.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 5.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 5.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 6.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 6.fasta.gz",
           "./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 7.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 7.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 7.fasta.gz")

for (file in seq_along(files)) {
  seg <- read.dna(gzfile(files[file]), "fasta")
  segPD <- phyDat(seg, type = "DNA", levels = NULL)
  mt <- modelTest(segPD)
  write.table(mt, file = paste(substring(files[file], 1, nchar(files[file]) - nchar(".fasta.gz")), "_mt.csv", sep=""), sep=",", quote = FALSE, row.names = F)
}

# Reconstruct maximum-likelihood trees under the chosen model of evolution.

files <- c("./Human influenza H3N2 segment 1-PB2 vRNA MSA - 2005-2014 OTU 1.fasta.gz", 
           "./Human influenza H3N2 segment 2-PB1 vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 3-PA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 4-HA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 5-NP vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 6-NA vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 7-M vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 segment 8-NS vRNA MSA - 2005-2014 OTU 1.fasta.gz",
           "./Human influenza H3N2 concatenated full-length vRNA MSA - 2005-2014 OTU 1.fasta.gz")

size <- 1000 # number of bootstrap replicates

for (file in seq_along(files)) {
  vRNA <- readDNAStringSet(files[file])
  d <- DistanceMatrix(vRNA, type="dist", correction="JC")
  tree <- IdClusters(d, method="ML", type="dendrogram", model="HKY85", myXStringSet=vRNA) # change 'model' parameter if desired.
  
  f <- function(x) {
    if (is.null(attributes(x)$leaf)) {   
      x0 <- paste(sort(unlist(x)), collapse=" ")   
      x1 <- f(x[[1]])                                 
      x2 <- f(x[[2]])
      return(list(x0, x1, x2))
    } else {
      return(NULL)
    }
  }
  
  # Perform bootstrapping.
  
  pBar <- txtProgressBar(style=3)
  bootstraps <- list()
  l <- unique(width(vRNA))
  for (i in seq_len(size)) { 
    r <- sample(l, replace=TRUE)
    at <- IRanges(r, width=1)            
    vRNA2 <- extractAt(vRNA, at)
    vRNA2 <- lapply(vRNA2, unlist)
    vRNA2 <- DNAStringSet(vRNA2)
    
    d <- DistanceMatrix(vRNA2, type="dist", correction="JC", verbose=FALSE)
    temp <- IdClusters(d, method="ML", type="dendrogram", model="HKY85", myXStringSet=vRNA2, verbose=FALSE)
    bootstraps[[i]] <- unlist(f(temp))   
    setTxtProgressBar(pBar, i/size)
  }
  
  bootstraps <- table(unlist(bootstraps))     
  original <- unlist(f(tree))       
  hits <- bootstraps[original]     
  names(hits) <- original
  w <- which(is.na(hits))
  if (length(w) > 0)         
    hits[w] <- 0
  hits <- round(hits/size*100)
  
  f <- function(x) {
    if (is.null(attributes(x)$leaf)) {
      attr(x, "edgetext") <- as.character(hits[paste(sort(unlist(x)), collapse=" ")])
    }
    return(x)
  }
  d <- dendrapply(tree, f)                  
  attr(d, "edgetext") <- NULL                
  
  WriteDendrogram(d, file = paste(substring(files[file], 1,  nchar(files[file]) - nchar(".fasta.gz")), "_MLtree", sep=""))
  par(mai= c(1, 0.5, 0.25, 0.1)) # bottom, left, top, right
  plot(d, edgePar=list(t.cex=0.5), nodePar=list(lab.cex=0.7, pch=NA), edge.root=FALSE)
  
}

# QUANTIFICATION OF TREE DISTANCES

# Trees built using the DECIPHER package must be read in by the ape package for analysis. 
# This converts the trees from objects of type = 'dendrogram' to type = 'phylo', 
# which introduces syntax errors that must be corrected prior to analysis.

# Working with one set of replicates at a time:

otu <- read.csv("./H3N2 2005-2014 Strains Analyzed - ape.csv") # Strain names reformatted without spacing for ape package.

subtype <- "H3N2"
segments <- c("1-PB2", "2-PB1", "3-PA", "4-HA", 
              "5-NP", "6-NA", "7-M", "8-NS")
nSegments <- length(segments)

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 1_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]])) # Replace strain names with cluster IDs in tree so that replicates can be compared.
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_1[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_1[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_1[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_1[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_1[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_1[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_1[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_1[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_1[9])
      treetips[j] = 9
    else if (trees[[i]]$tip.label[j] == otu$id_1[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_1[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_1[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  par(mai = c(0.5, 0.5, 0.5, 0.5))
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

# Combine all eight trees from each vRNA segment into one object and calculate tree distances.

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf1 <- dist.topo(unroot(trees), method = "PH85") # Robinson-Foulds distance
cid1 <- TreeDistance(trees) # Clustering Information Distance

# Repeat for remaining replicates.

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 2_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_2[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_2[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_2[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_2[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_2[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_2[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_2[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_2[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_2[9])
      treetips[j] = 9
    else if (trees[[i]]$tip.label[j] == otu$id_2[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_2[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_2[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf2 <- dist.topo(unroot(trees), method = "PH85")
cid2 <- TreeDistance(trees)

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 3_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_3[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_3[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_3[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_3[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_3[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_3[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_3[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_3[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_3[9])
      treetips[j] = 9    
    else if (trees[[i]]$tip.label[j] == otu$id_3[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_3[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_3[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf3 <- dist.topo(unroot(trees), method = "PH85")
cid3 <- TreeDistance(trees)

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 4_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]])) 
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_4[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_4[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_4[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_4[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_4[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_4[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_4[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_4[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_4[9])
      treetips[j] = 9
    else if (trees[[i]]$tip.label[j] == otu$id_4[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_4[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_4[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf4 <- dist.topo(unroot(trees), method = "PH85")
cid4 <- TreeDistance(trees)

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 5_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_5[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_5[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_5[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_5[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_5[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_5[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_5[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_5[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_5[9])
      treetips[j] = 9
    else if (trees[[i]]$tip.label[j] == otu$id_5[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_5[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_5[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf5 <- dist.topo(unroot(trees), method = "PH85")
cid5 <- TreeDistance(trees)

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 6_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_6[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_6[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_6[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_6[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_6[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_6[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_6[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_6[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_6[9])
      treetips[j] = 9
    else if (trees[[i]]$tip.label[j] == otu$id_6[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_6[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_6[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  plot(trees[[i]])
  add.scale.bar(0, 0.77, length = 0.005, col = "black")
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="red")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf6 <- dist.topo(unroot(trees), method = "PH85")
cid6 <- TreeDistance(trees)

files <- paste("./Human influenza", subtype, "segment", segments, "vRNA MSA - 2005-2014 OTU 7_MLtree")

length <- length(files)
trees <- vector(mode = "list", length = length)

for (i in seq_along(trees)) {
  
  trees[[i]] <- read.tree(files[i])
  trees[[i]]$node.label <- gsub("\"", "", trees[[i]]$node.label)
  trees[[i]]$tip.label <- gsub("\"", "", trees[[i]]$tip.label)
  
  treetips <- integer(Ntip(trees[[i]]))
  for (j in seq_along(treetips)) {
    if (trees[[i]]$tip.label[j] == otu$id_7[1])
      treetips[j] = 1
    else if (trees[[i]]$tip.label[j] == otu$id_7[2])
      treetips[j] = 2
    else if (trees[[i]]$tip.label[j] == otu$id_7[3])
      treetips[j] = 3
    else if (trees[[i]]$tip.label[j] == otu$id_7[4])
      treetips[j] = 4
    else if (trees[[i]]$tip.label[j] == otu$id_7[5])
      treetips[j] = 5
    else if (trees[[i]]$tip.label[j] == otu$id_7[6])
      treetips[j] = 6
    else if (trees[[i]]$tip.label[j] == otu$id_7[7])
      treetips[j] = 7
    else if (trees[[i]]$tip.label[j] == otu$id_7[8])
      treetips[j] = 8
    else if (trees[[i]]$tip.label[j] == otu$id_7[9])
      treetips[j] = 9
    else if (trees[[i]]$tip.label[j] == otu$id_7[10])
      treetips[j] = 10
    else if (trees[[i]]$tip.label[j] == otu$id_7[11])
      treetips[j] = 11
    else if (trees[[i]]$tip.label[j] == otu$id_7[12])
      treetips[j] = 12
    else (treetips[j] = NA)
  }
  
  trees[[i]]$tip.label <- treetips
  plot(trees[[i]], font = 1, edge.width = 3, cex = 1.1)
  add.scale.bar(cex = 1.2, length = 0.01, ask = T, font = 2)
  nodelabels(trees[[i]]$node.label, adj = c(1.2, -0.5), frame = "none", cex=0.8, col="darkred")
}

trees <- c(trees[[1]], trees[[2]], trees[[3]], trees[[4]], trees[[5]], trees[[6]], trees[[7]], trees[[8]])
rf7 <- dist.topo(unroot(trees), method = "PH85")
cid7 <- TreeDistance(trees)

lr <- list(as.matrix(rf1), as.matrix(rf2), as.matrix(rf3), as.matrix(rf4), as.matrix(rf5), as.matrix(rf6), as.matrix(rf7))
rf <- do.call(cbind, lr)
rf <- array(rf, dim=c(dim(lr[[1]]), length(lr)))

lc <- list(as.matrix(cid1), as.matrix(cid2), as.matrix(cid3), as.matrix(cid4), as.matrix(cid5), as.matrix(cid6), as.matrix(cid7))
cid <- do.call(cbind, lc)
cid <- array(cid, dim=c(dim(lc[[1]]), length(lc)))

r2 <- apply(rf, c(1, 2), mean) # Compute the mean Robinson-Foulds distances across replicates.
row.names(r2) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
colnames(r2) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

rsem2 <- apply(rf, c(1, 2), std.error) # Compute the standard error.
write.table(r2, file = "./H3N2 2005-2014 ML Tree  RF_mean.csv", sep=",", quote = FALSE, row.names = T)
write.table(rsem2, file = "./H3N2 2005-2014 ML Tree  RF_sem.csv", sep=",", quote = FALSE, row.names = T)

c2 <- apply(cid, c(1, 2), mean) # Compute the mean Robinson-Foulds distances across replicates.
row.names(c2) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
colnames(c2) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

csem2 <- apply(cid, c(1, 2), std.error) # Compute the standard error.
write.table(c2, file = "./H3N2 2005-2014 ML Tree  CID_mean.csv", sep=",", quote = FALSE, row.names = T)
write.table(csem2, file = "./H3N2 2005-2014 ML Tree  CID_sem.csv", sep=",", quote = FALSE, row.names = T)

# The data can be visualized in a number of ways. Here is how to generate a distance network.

network <- matrix(r2, nrow = 8)
row.names(network) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
colnames(network) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
tree <- IdClusters(network, method="UPGMA", type="dendrogram")

par(mai = c(1, 1, 0.75, 0.25))  # bottom, left, top, right
plot(tree, xlab = "vRNA Segment", ylab = "Robinson-Foulds distance", type = "triangle", 
     cex.lab = 1.2, cex= 1.2, main = "2005-2014", ylim= c(0, 8)) 
WriteDendrogram(tree, file = "./H3N2 2005-2014 RF network") 

network <- matrix(c2, nrow = 8)
row.names(network) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
colnames(network) <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")
tree <- IdClusters(network, method="UPGMA", type="dendrogram")

plot(tree, xlab = "vRNA Segment", ylab = "Clustering Information Distance", type = "triangle", 
     cex.lab = 1.2, cex= 1.2, main = "2005-2014", ylim= c(0, 0.25)) 
WriteDendrogram(tree, file = "./H3N2 2005-2014 CID network") 

# Here is how to generate a distance heatmap.

# First, remove the duplicated portions of the data.

r2 <- r2[,-8] # Remove the NS column.
r2 <- r2[-1,] # Remove the PB2 row.
r2[1,][2:7] <- NA
r2[2,][3:7] <- NA
r2[3,][4:7] <- NA
r2[4,][5:7] <- NA
r2[5,][6:7] <- NA
r2[6,][7] <- NA

c2 <- c2[,-8] # Remove the NS column.
c2 <- c2[-1,] # Remove the PB2 row.
c2[1,][2:7] <- NA
c2[2,][3:7] <- NA
c2[3,][4:7] <- NA
c2[4,][5:7] <- NA
c2[5,][6:7] <- NA
c2[6,][7] <- NA

rsem2 <- rsem2[,-8]
rsem2 <- rsem2[-1,]
rsem2[1,][2:7] <- NA
rsem2[2,][3:7] <- NA
rsem2[3,][4:7] <- NA
rsem2[4,][5:7] <- NA
rsem2[5,][6:7] <- NA
rsem2[6,][7] <- NA

csem2 <- csem2[,-8]
csem2 <- csem2[-1,]
csem2[1,][2:7] <- NA
csem2[2,][3:7] <- NA
csem2[3,][4:7] <- NA
csem2[4,][5:7] <- NA
csem2[5,][6:7] <- NA
csem2[6,][7] <- NA

rownames(rsem2) <- as.character(c("PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(rsem2) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M"))

rownames(csem2) <- as.character(c("PB1", "PA", "HA", "NP", "NA", "M", "NS"))
colnames(csem2) <- as.character(c("PB2", "PB1", "PA", "HA", "NP", "NA", "M"))

palette <- colorRampPalette(brewer.pal(7, "YlGnBu"))(20)
cols <- palette[20:1]

par(mai = c(0.5, 0.5, 1, 1))
levelplot(r2, 
          at=c(0, 2, 4, 6, 8, 10, 12, 14, 16),
          xlab = "",
          ylab = "",
          col.regions = cols)

levelplot(rsem2, 
          at=c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
          xlab = "",
          ylab = "",
          col.regions = cols)

levelplot(c2, 
          at=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6),
          xlab = "",
          ylab = "",
          col.regions = cols)

levelplot(csem2, 
          at=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,0.07, 0.08, 0.09, 0.1),
          xlab = "",
          ylab = "",
          col.regions = cols)













