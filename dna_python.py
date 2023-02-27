from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import click


# GC content function

def gc_content(sequence):
    print()
    print('The GC content is..')

    
    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()                                     # Remove trailing white space

    # Calculating GC content
    bac=round(gc_fraction(seq)*100,2)                           # gc_fraction is function imported to calculate GC content.
    print('GC ',bac,'%')


# Dinucleotide content

def di_nuc(sequence):
    print()
    print('Dinucleotide frequencies..')
    
    
    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()

    AA=round(100*seq.count_overlap('AA')/len(seq),2)            # count_overlap method counts AA in a sequence. overlap method count 3 for AA in 'AAA', whereas only .count method counts 1 AA for 'AAA'.
    TT=round(100*seq.count_overlap('TT')/len(seq),2)
    GG=round(100*seq.count_overlap('GG')/len(seq),2)
    CC=round(100*seq.count_overlap('CC')/len(seq),2)


    print('AA',AA,"%")
    print('TT',TT,"%")
    print('GG',GG,"%")
    print('CC',CC,"%")


# Trinucleotide content

def tri_nuc(sequence):
    print()
    print('Trinucleotide frequencies..')

    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()

    AAA=round(100*seq.count_overlap('AAA')/len(seq),2)            # count_overlap method counts AA in a sequence. overlap method count 3 for AA in 'AAA', whereas only .count method counts 1 AA for 'AAA'.
    TTT=round(100*seq.count_overlap('TTT')/len(seq),2)
    GGG=round(100*seq.count_overlap('GGG')/len(seq),2)
    CCC=round(100*seq.count_overlap('CCC')/len(seq),2)


    print('AAA',AAA,"%")
    print('TTT',TTT,"%")
    print('GGG',GGG,"%")
    print('CCC',CCC,"%")


# Ambiguous Frequency

def ambi(sequence):
    print()
    print('Ambiguous frequencies..')

    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()

    NN=round(100*seq.count_overlap('N')/len(seq),2)

    print('N',NN,"%")



# Command Line Arguments
# click 

@click.command()                                                            # Using click to create options.
@click.option('-g', is_flag=True,help='Calculate GC content')	            # is_flag is set to true, so when -g is mentioned, the code under it will run.
@click.option('-d',is_flag=True,help='Calculate dinucleotide frequency')
@click.option('-t',is_flag=True,help='Calculate trinucleotide frequency')
@click.option('-n',is_flag=True,help='Calculate ambiguous frequency')
@click.argument('filename',type=click.Path(exists=True))                    # click.path adds functionality to check if the filepath exists. 
def data(g,d,t,n,filename):
    # Adding a docstring
    """
    Calculate GC content, Dinucleotide, Trinucleotide and Ambiguous Frequencies.

    Input Arguments:
    Filename

    """
    # Reading Fasta File
    sequence = [i for i in SeqIO.parse(filename,'fasta')]      # Parsing fasta file and storing each sequence inside a list

    if g:
        gc_content(sequence)
    if d:
        di_nuc(sequence)
    if t:                           # elif is not used, because once elif is true python stops interpreting further conditions
        tri_nuc(sequence)
    if n:
        ambi(sequence)
    

if __name__=='__main__':            # execute data() when the file is run directly and is not imported
    data()