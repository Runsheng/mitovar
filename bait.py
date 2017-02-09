#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/8 17:00
# @Author  : Runsheng     
# @File    : bait.py

"""
functions used to bait the reads from mitochondrial out
"""

from utils import myexe
import os
import pysam
from glob import glob


def bwa_index_wrapper(ref_file):
    """
    :param ref_file:
    :return:
    """
    cmd_index="bwa index {ref}".format(ref=ref_file)
    myexe(cmd_index)
    return ref_file


def sra2fq_wrapper(sra_file, outdir):

    sra_cmd = "fastq-dump --split-files {} --outdir {}".format(sra_file, outdir)
    print(sra_cmd)
    print(myexe(sra_cmd))

    return outdir


def get_fq_dict(workdir_fastq):
    """
    :param workdir_fastq:
    :return: a dict with {"ERRxxxx": ["ERRxxxx_1.fq","ERRxxxx_2.fq"], "ERRxxx2":["ERRxxx2_1.fq"]}
    """
    fq_dict={}
    os.chdir(workdir_fastq)
    fq_names=glob("*.fq*")+glob("*.fastq*")
    for name in fq_names:
        if ".fq" or ".fastq" in name:
            prefix=name.split("_")[0]
            name_abs=os.path.join(workdir_fastq, name)
            try:
                fq_dict[prefix].append(name_abs)
            except KeyError:
                fq_dict[prefix]=[]
                fq_dict[prefix].append(name_abs)
    return fq_dict


def bwa_mem_wrapper(ref_file, fq_str, core=25, min_seed_length=20, band_width=2000, out="mapped.bam"):
    """
    a bwa mem mapper only collect the mapped reads
    :param ref_file:
    :param fq_str: a list contains the name of the sra -extracted
    :param core, min_seed_length and band_width is bwa mem parameter -t, -k and -w, respectively
    :return: the out bam file
    """

    cmd_bwa="bwa mem -k {band} -w {width} -t {core} {ref} \
        {fq_str} \
        | samtools view -F 4  -b -o {out}".format(
        band=min_seed_length, width=band_width, core=core,ref=ref_file,
        fq_str=fq_str,
        out=out)
    print(cmd_bwa)
    myexe(cmd_bwa)
    return out


def sort_index_wrapper(bamfile, core=1, bam_sorted=None):
    if bam_sorted==None:
        bam_sorted=bamfile.split(".")[0]+"_s.bam"
    else:
        pass

    sort_cmd="samtools sort {bamfile} -@ {core} -o {bam_sorted}".format(
        bamfile=bamfile, core=core, bam_sorted=bam_sorted)
    index_cmd="samtools index {bam_sorted}".format(bam_sorted=bam_sorted)

    myexe(sort_cmd)
    myexe(index_cmd)
    return bam_sorted


def get_name(bam_sorted):
    """
    from the bam file get the mapped read name
    :param bam_sorted:
    :return:
    """

    samfile=pysam.AlignmentFile(bam_sorted, "rb")
    names_set=set()
    for record in samfile:
        names_set.add(record.qname)

    return names_set


def write_name(names_set, namefile="name.lst"):
    """
    used in fq_subset_main(
    write a namefile and return the name of this file
    """
    with open(namefile, "w") as fw:
        for name in names_set:
            fw.write(name)
            fw.write("\n")
    return namefile


def fq_subseq(fqin, namefile, fqout=None):
    """
    used in fq_subset_main()
    :param fqin:454
    :param namefile:
    :param fqout:
    :return:
    """
    if fqout==None:
        fqout=fqin.split(".")[0]+"_s.fq"

    cmd_subseq = "seqtk subseq {fqin} {namefile} > {fq_out}".format(fqin=fqin, namefile=namefile, fq_out=fqout)
    print(cmd_subseq)
    myexe(cmd_subseq)
    return fqout


def fq_subset_main(fqin, names_set, fqout=None):
    namefile = write_name(names_set, namefile="name.lst")
    fqout=fq_subseq(fqin, namefile, fqout)
    return fqout



