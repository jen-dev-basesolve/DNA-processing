#########################################
# GC content function
#########################################

def gc_content(file_name):
    import Bio

    # Reading Fasta File
    from Bio import SeqIO
    sequence = [i for i in SeqIO.parse(file_name,'fasta')]
    seq=''
    for i in sequence:
        seq+=i.seq.rstrip()

    # Calculating GC content
    from Bio.SeqUtils import gc_fraction
    bac=round(gc_fraction(seq)*100,2)
    print(bac,'%')


############################################
# Command Line Arguments
# click 
##############################

import click

@click.command()
@click.option('-g', is_flag=True,help='Type -g to calculate GC content')		# is_flag is set to true, so when -g is mentioned, the code under it will run.
@click.option('-d',is_flag=True)
@click.option('-t',is_flag=True)
@click.option('-n',is_flag=True)
def data(g,d,t,n):
    if g:
        gc_content()
    if d:
        print('This is d')
    if t:                               # elif is not used, because once elif is true, then python stops interpreting
        print('This is t')
    if n:
        print('This is n')
    else:
        print('error')

if __name__=='__main__':
    data()