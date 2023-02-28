# GC content and Dinucleotide frequency calculations

library(Biostrings)
library(LncFinder)
library(seqinr)

filename <- "MonkeyPox.fasta"

# GC content
gc_count <- function(filename) 
{
  samfa <- read.fasta(file=filename)	# File should be in current working directory
  compute_GC(samfa,parallel.cores = 2)		# parallel.cores=2 uses two cores for computation, takes less time.
}

# Dinucleotide frequency

dinuc_freq <- function(filename)
{samfa <- read.fasta(file=filename)	# File should be in current working directory

for (i in 1:length(samfa))
{
  dinuc_table <- count(samfa[[i]],freq=TRUE,wordsize=2)
  aa <- dinuc_table[1]
  tt <- dinuc_table[16]
  gg <- dinuc_table[11]
  cc <- dinuc_table[6]
  cat("SeqNo.",i,'\n')
  cat('aa=',aa,"tt=", tt,"cc=", cc, "gg=", gg,'\n' )
}

}

gc_count(filename)
dinuc_freq(filename)

###############################
# Issues: - Does not give cumulative output for all sequences
#         - Does not display sequence ID before each output for Dimer frequency