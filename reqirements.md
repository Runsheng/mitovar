### The third-part software used
#### Bins need to be in the $PATH
1. BWA-MEM
2. Spades assembler >=3.9
3. seqtk
4. exonerate: for mito gene annotation, included in ./bins
5. MiTFi package: for tRNA and rRNA annotation, included in ./bins/mitfi
6. infernal=1.0.2, required by MiTFi


#### python packages
1. Biopython
2. pysam

#### install dependencies 
The easiest way to install all dependencies is conda throught bioconda.
```
conda install -c bioconda \
    biopython pysam seqtk bwa infernal=1.0.2 spades 
```

#### optional packages
1. fastq-dump (from sra-toolkit): if start from sra file other then fastq