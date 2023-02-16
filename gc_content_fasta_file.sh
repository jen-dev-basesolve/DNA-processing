#! /bin/bash
filename=$@                     # collecting all the filenames as parameters
for f in $filename              # Looping over files    
do
    echo " $f is being processed..."                            
    gc=( $( grep -v ">" < "$f" | grep -io 'g\|c'< "$f" | wc -l))    # Reading lines that dont start with < using -v. grep -io matches to either g or c and outputs each match on single line. wc -l counts the number of lines or indirectly the number of g and c. This is stored in a variable. 
    total=( $( grep -v ">" < "$f" | tr -d '\s\r' | wc -c))          # Spaces, tabs, new line are removed from the file using tr. Then the number of characters are counted by wc -c
    percent=( $( echo "scale=2;100*$gc/$total" |bc -l))             # bc -l is used to get the answer in float format. scale=2 mentions the number of decimal points.
    echo " The GC content of $f is: "$percent"%"
    echo
done 