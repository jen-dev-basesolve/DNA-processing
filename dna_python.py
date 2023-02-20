# Read the file from user input
file_name=input('Enter File name:')  # File input can further be improved using getopt.getopt()
print('Loading...')

#########################################
# Help
#########################################

def help():
    print ('Python script to calculate GC content')
    print (' g - GC content \n h - Help')



#########################################
# GC content function
#########################################

def gc_content(file_name):
    # Reading lines that don't start with '>'

    data=''                             # Creating and empty string
    with open(file_name,'r') as file:
        for line in file.readlines():
            if not line.startswith('>'):
                data=data+line.rstrip() # rstrip removes any whitespace

    # Calculating number of G and C

    gc_count=0
    for i in data:
        if i == 'C' or i == 'G':
            gc_count=gc_count+1

    print(f'The total count of GC content is: {gc_count}')
    
    # Total number of bases
    seq_len=len(data)       

    # GC percent
    gc_per=round(100*(gc_count/seq_len),2)
    print("The percentage of GC content is:",gc_per,"%")



############################################
# Command Line Arguments
############################################

import sys

for i in sys.argv:
    if i=='g':
        gc_content(file_name)
    elif i=='h': 
        help()
