# Mitovar
### Mitochondrial DNA sequence assembly and annotation for nematodes using NGS data. 
- Dataset S3: in folder ./data 
- Stable release: V0.99 [link]



### Install
##### Install dependencies 
The easiest way to install all dependencies is conda throught bioconda.
```bash
conda install -c bioconda \
    biopython pysam seqtk bwa infernal=1.0.2 spades 
```

##### Install mitovar (with mitfi and exonerate)
```bash
git clone https://github.com/Runsheng/mitovar.git
```

##### Add mitovar.py to $PATH
```bash
export PATH="mitovar_install_path":$PATH
```


##### Optional 
- fastq-dump (from sra-toolkit): if start from sra file other than fastq
- The pre-compiled exonerate only works on Linux x64 system, tested in Ubuntu 16.04, 14.04 and CentOS 7.

### Run commands

#### The command contains:
- anno: annotate a mtDNA fasta sequences and generate a tbl file for genbank submission
- assemble: assemble the mtDNA fasta file from a NGS fastq file and a nearby reference file


An example to run annotation command:
```bash
mitovar.py anno -f mtDNA.fasta -c cel_p.fa -r cel_rrna.fa -s cel
```

An example to run assemble command:

*Note: put all fastq files in ./{spe}/fastq, and run command in ./*

*For instance, the species name {spe} is cbr*

*run the following command in /home/user, the fastq should be in /home/usr/cbr/fastq*
```bash
mitovar.py assemble -f cel.fa -p 32 -s cbr
```