#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/8 17:00
# @Author  : Runsheng     
# @File    : bait.py

from utils import myexe
import os
import pysam

def bwa_index_wrapper(ref_file):
    """
    :param ref_file:
    :return:
    """
    cmd_index="bwa index {ref}".format(ref_file)
    return ref_file


def sra2fq_wrapper(srafile, outdir):

    sra_cmd = "fastq-dump --split-files {} --outdir {}".format(srafile, outdir)

    os.list(outdir)


def bwa_mem_wrapper(ref_file, fastq_f, fastq_r, core=25, min_seed_length=18, band_width=2000, out="mapped.bam"):
    """
    a bwa mem mapper only collect the mapped reads
    :param ref_file:
    :param fastq_f:
    :param fastq_r:
    :param core, min_seed_length and band_width is bwa mem parameter -t, -k and -w, respectively
    :return: the out bam file
    """

    cmd_bwa="bwa mem -k {band} --w {width} -t {core} {ref} \
        {fastq_f} {fastq_r} \
        | samtools view -F 4  -b -o {out}".format(
        band=min_seed_length, width=band_width, core=core,ref=ref_file,
        fastq_f=fastq_f, fastq_r=fastq_r,
        out=out)
    print(cmd_bwa)
    myexe(cmd_bwa)
    return out


def sort_index_wrapper(bamfile, core=1, bam_sorted=None):
    if bam_sorted==None:
        bam_sorted=bamfile.split(".")[0]+"_s.bam"
    else:
        pass

    sort_cmd="samtools sort bamfile -@ {core} -o {bam_sorted}".format(
        core=core, bam_sorted=bam_sorted)
    index_cmd="samtools index {bam_sorted}".format(bam_sorted=bam_sorted)

    myexe(sort_cmd)
    myexe(index_cmd)
    return bam_sorted


def get_name(bam_sorted):

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



