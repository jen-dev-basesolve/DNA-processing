# import libraries

library(Biostrings)
library(LncFinder)
library(seqinr)
library(optparse)


# GC content function
gc_count <- function(filename) 
{
  # read.fasta{seqinr} : pares fasta file into biostrings, similar to lists
  samfa <- read.fasta(file=filename)	# File should be in current working directory
  # compute_GC{LncFinder} : Compute GC count, with overlaping
  compute_GC(samfa)		# parallel.cores=2 uses two cores for computation, takes less time.
}


# Dinucleotide frequency function
dinuc_freq <- function(data,tot_nucs)
{
  # vcountPattern{Biostrings} : count pattern over multiple xstrings
  aa <- vcountPattern(subject=data,pattern="AA",max.mismatch = 0)
  aa <- 100*sum(aa)/tot_nucs
  
  tt <- vcountPattern(subject=data,pattern="TT",max.mismatch = 0)
  tt <- 100*sum(tt)/tot_nucs
  
  gg <- vcountPattern(subject=data,pattern="GG",max.mismatch = 0)
  gg <- 100*sum(gg)/tot_nucs
  
  cc <- vcountPattern(subject=data,pattern="CC",max.mismatch = 0)
  cc <- 100*sum(cc)/tot_nucs
  
  
  cat("AA",aa,"% \n")
  cat("TT",tt,"% \n")
  cat("GG",gg,"% \n")
  cat("CC",cc,"% \n")
}

# Trinucleotide frequency
trinuc_freq <- function(data,tot_nucs)
{
  
  aa <- vcountPattern(subject=data,pattern="AAA",max.mismatch = 0)
  aa <- 100*sum(aa)/tot_nucs
  
  tt <- vcountPattern(subject=data,pattern="TTT",max.mismatch = 0)
  tt <- 100*sum(tt)/tot_nucs
  
  gg <- vcountPattern(subject=data,pattern="GGG",max.mismatch = 0)
  gg <- 100*sum(gg)/tot_nucs
  
  cc <- vcountPattern(subject=data,pattern="CCC",max.mismatch = 0)
  cc <- 100*sum(cc)/tot_nucs
  
  cat("\n")
  cat("AAA",aa,"% \n")
  cat("TTT",tt,"% \n")
  cat("GGG",gg,"% \n")
  cat("CCC",cc,"% \n")
}

# Ambiguous frequency
ambi_freq <- function(data,tot_nucs)
{
  nn <- vcountPattern(subject=data,pattern="N",max.mismatch = 0)
  nn <- 100*sum(nn)/tot_nucs
  cat("\n")
  cat("N",nn,"% \n")
}


# Main Program
main <- function(){
  # make_option{Optparse} : Similar to python tool used to add options to Rscripts
  option_list = list(
    # Takes long form and short form options.
    # action="store" instructs to store user input after option
    make_option(c("-f", "--file"), action="store", default=NA, type='character',
                help="Enter filename or path"),
    # action="store_true" : store True when the option is called in cli
    make_option(c("-g", "--gc_content"), action="store_true", default=FALSE,
                help="Calculate GC content"),  
    make_option(c("-d", "--di_nuc"), action="store_true", default=FALSE,
                help="Calculate Dinucleotide Frequencies"),
    make_option(c("-t", "--tri_nuc"), action="store_true", default=FALSE,
                help="Calculate Trinucleotide Frequencies"),
    make_option(c("-a","--ambiguous"), action="store_true", default=FALSE,
                help="Calculate Ambiguous nucleotide Frequencies")
  )
  # parse_args : A method to analysize cli option and take appropiate action
  
  opt <- parse_args(OptionParser(option_list=option_list))
  
  # parsing fasta file to Xstrings 
  data <- readDNAStringSet(filepath = opt$f, format="fasta")

  # calculate total number of nucleotides in fasta file
  tot_nucs=0
  for ( i in 1:length(data))
  {
    tot_nucs=tot_nucs+length(data[[i]])
  }
  
  # command line conditions

  if (is.null(opt$f)){
    print("ERROR, please enter Filename")
  }
  
  if (opt$g){
    gc_count(opt$f)
  }
  
  if (opt$d){
    dinuc_freq(data,tot_nucs)
  }
  
  if (opt$t){
    trinuc_freq(data,tot_nucs)
  }
  
  if (opt$a){
    ambi_freq(data,tot_nucs)
  }
  
}

main()