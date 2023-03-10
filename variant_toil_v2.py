import os

# Downloading Reference Genome of ecoli
def ref_genome():
    comms = "mkdir dc_workshop_1 ; cd dc_workshop_1 ; mkdir -p data/ref_genome "
    os.system(comms)
    comms = "curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
    os.system(comms)
    comms = "gunzip data/ref_genome/ecoli_rel606.fasta.gz"
    os.system(comms)


# Downloading trimmed FastQ files for faster operations
def fastq_files():
    comms = "cd dc_workshop_1 ; curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248 ; tar xvf sub.tar.gz ; mkdir data/trimmed_fastq_small ; mv sub/ data/trimmed_fastq_small"
    os.system(comms)
    comms = "mkdir -p results/sam results/bam results/bcf results/vcf"
    os.system(comms)


# Indexing reference fasta file
def index_ref():
    comms = "singularity exec bwa_latest.sif bwa index dc_workshop_1/data/ref_genome/ecoli_rel606.fasta"
    os.system(comms)


# Aligning sample sequences to reference genome
def align_fastq():
    comms = ("singularity exec bwa_latest.sif bwa mem dc_workshop_1/data/ref_genome/ecoli_rel606.fasta dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_1.trim.sub.fastq dc_workshop_1/data/trimmed_fastq_small/sub/SRR2584866_2.trim.sub.fastq > dc_workshop_1/results/sam/SRR2584866.aligned.sam")
    os.system(comms)


# Converting SAM file to BAM file using view option in samtools
# Sorting bam files
def convrt_sort():
    comms = "singularity exec samtools_latest.sif samtools view -S -b dc_workshop_1/results/sam/SRR2584866.aligned.sam > dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    os.system(comms)
    comms = "singularity exec samtools_latest.sif samtools sort -o dc_workshop_1/results/bam/SRR2584866.aligned.sorted.bam dc_workshop_1/results/bam/SRR2584866.aligned.bam"
    os.system(comms)
