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
```
export PATH="mitovar_install_path":$PATH
```


##### Optional 
fastq-dump (from sra-toolkit): if start from sra file other than fastq


### Run commands

#### The command contains:
- anno: annotate a mtDNA fasta sequences and generate a tbl file for genbank submission
- assemble: assemble the mtDNA fasta file from a NGS fastq file and a nearby reference file


A example for running annotation command:
```
mitovar.py anno -f mtDNA.fasta -c cel_p.fa -r cel_rrna.fa -s cel
mitovar.py assemble -f cel.fa -p 32 -s cbr
```