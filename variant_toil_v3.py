from toil.common import Toil
from toil.job import Job
import os
import tempfile
import requests
import shlex
import subprocess


def helloWorld(message,memory="2G", cores=2, disk="3G"):
    return message

# def ref_genome(memory="2G", cores=2, disk="3G"):
#     parent_dir = "/home/bioinfo/singularity/variant_analysis"
#     child_dir = "dc_workshop_1/data/ref_genome"
#     path = os.path.join(parent_dir,child_dir)
#     os.makedirs(path,exist_ok=True)
#     # Downloading fasta file
#     url = "http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"
#     local_filename = url.split('/')[-1]
#     # NOTE the stream=True parameter below
#     with requests.get(url, stream=True) as r:
#         r.raise_for_status()
#         with open("/home/bioinfo/singularity/variant_analysis/dc_workshop_1/data/ref_genome/%s"%local_filename, 'wb') as f:
#             for chunk in r.iter_content(chunk_size=8192): 
#                 f.write(chunk)

    # # extracting fasta file
    # comms = ("gunzip %s/%s"%(path,local_filename))
    # os.system(comms)    


# Downloading trimmed FastQ files for faster operations
def fastq_files():
    parent_dir = "/home/bioinfo/singularity/variant_analysis"
    child_dir = "dc_workshop_1/data/trimmed_fastq_small"
    path = os.path.join(parent_dir,child_dir)
    os.makedirs(path,exist_ok=True)
    
    comms = ("curl -L -o %s/sub.tar.gz https://ndownloader.figshare.com/files/14418248"%path)
    args = shlex.split(comms)
    print(subprocess.run(args,shell=False,capture_output=True))  # run correct
    
    comm2 = ("tar xvf /home/bioinfo/singularity/variant_analysis/dc_workshop_1/data/trimmed_fastq_small/sub.tar.gz")
    os.system(comm2)
    # args = shlex.split(comms)
    # print(subprocess.run(args,shell=False,capture_output=True))

    # comms = "cd %s/dc_workshop_1 ; mkdir -p results/sam results/bam results/bcf results/vcf"%parent_dir
    # os.system(comms)
    

if __name__=='__main__':
    jobstore: str = tempfile.mkdtemp("test")
    os.rmdir(jobstore)
    options = Job.Runner.getDefaultOptions(jobstore)
    options.logLevel='DEBUG'
    j1 = Job.wrapFn(helloWorld,message='Hi There >>>>>>>>>>')
#    j2 = Job.wrapFn(ref_genome)
    j3 = Job.wrapFn(fastq_files)

    with Toil(options) as toil:
        print(toil.start(j1))
#        print(toil.start(j2))
        print(toil.start(j3))