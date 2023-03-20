#! /bin/bash

######################################
# Help
######################################
function Help()
{
    # Display help
    echo
    echo "################## HELP #########################"
    echo
    echo "This script can be used to find gc content, dinucleotide percent, trinucleotide percent, ambigious percent"
    echo
    echo "Syntax:"
    echo "file.sh -f <file.fasta> -g|d||t|n"
    echo
    echo "Options:"
    echo "f- Input Fasta File"
    echo "g- GC content"
    echo "d- Dinucleotide percent"
    echo "t- Dinucleotide percent"
    echo "n- Ambigious N base pair percent"
    echo
    echo
}

########################################################
##################### Main Program #####################
########################################################

######################################
# GC Content Calculation
######################################

function gc_content()
{
echo
echo " $f is being processed..."                            
gc=( $( grep -v ">" < "$f" | grep -io 'g\|c'< "$f" | wc -l))    # Reading lines that dont start with < using -v. grep -io matches to either g or c and outputs each match on single line. wc -l counts the number of lines or indirectly the number of g and c. This is stored in a variable. 
total=( $( grep -v ">" < "$f" | tr -d '\s\r' | wc -c))          # Spaces, tabs, new line are removed from the file using tr. Then the number of characters are counted by wc -c
percent=( $( echo "scale=2;100*$gc/$total" |bc -l))             # bc -l is used to get the answer in float format. scale=2 mentions the number of decimal points.
echo " The GC content of $f is: "$percent"%"
echo
}

######################################
# Dinucleotide Frequency Calculation
######################################

function di_nuc()
{
echo " Dinucleotide frequencies are as follows.."

AA=( $( grep -v '>' < "$f" | grep -io 'AA'< "$f" | wc -l))
TT=( $( grep -v '>' < "$f" | grep -io 'TT'< "$f" | wc -l))
CC=( $( grep -v '>' < "$f" | grep -io 'CC'< "$f" | wc -l))
GG=( $( grep -v '>' < "$f" | grep -io 'GG'< "$f" | wc -l))

Total=( $( grep -v '>' < "$f" | tr -d '\n\r' < "$f" | wc -c))

PA=( $(echo "scale=2;100*$AA/$Total" | bc -l))
PT=( $(echo "scale=2;100*$TT/$Total" | bc -l))
PC=( $(echo "scale=2;100*$CC/$Total" | bc -l))
PG=( $(echo "scale=2;100*$GG/$Total" | bc -l))

echo " AA  $PA%" 
echo " TT  $PT%" 
echo " CC  $PC%" 
echo " GG  $PG%" 
}

######################################
# Trinucleotide Frequency Calculation
######################################

function tri_nuc()
{
echo
echo " Trinucleotide frequencies are as follows.."
echo

AAA=( $( grep -v '>' < "$f" | grep -io 'AAA'< "$f" | wc -l))
TTT=( $( grep -v '>' < "$f" | grep -io 'TTT'< "$f" | wc -l))
CCC=( $( grep -v '>' < "$f" | grep -io 'CCC'< "$f" | wc -l))
GGG=( $( grep -v '>' < "$f" | grep -io 'GGG'< "$f" | wc -l))

Total=( $( grep -v '>' < "$f" | tr -d '\n\r' < "$f" | wc -c))

PAA=( $(echo "scale=2;100*$AAA/$Total" | bc -l))
PTT=( $(echo "scale=2;100*$TTT/$Total" | bc -l))
PCC=( $(echo "scale=2;100*$CCC/$Total" | bc -l))        # Why is this part now showing 0 ahead of decimal?
PGG=( $(echo "scale=2;100*$GGG/$Total" | bc -l))

echo " AAA  $PAA%" 
echo " TTT  $PTT%" 
echo " CCC  $PCC%" 
echo " GGG  $PGG%" 
echo
}

#####################################
# Ambiguous Frequency
#####################################
function n_freq()
{
echo    
echo " Ambiguous frequency is..."
N=( $( grep -v '>' < "$f" | grep -io 'N'< "$f" | wc -l))
Total=( $( grep -v '>' < "$f" | tr -d '\n\r' < "$f" | wc -c))
PN=( $(echo "scale=2;100*$N/$Total" | bc -l))
echo " N  $PN%"
}

#####################################
# GET OPTIONS
#####################################

while getopts ntdhgf: option; do
	case $option in
	f) f=$OPTARG;;
    g) gc_content;;
    d) di_nuc;;
    t) tri_nuc;;
    n) n_freq;;
    h) Help; exit;;
    \?) echo "Invalid Option" ; Help exit ;;					# '*' also works in place of '\?'
	esac
done
