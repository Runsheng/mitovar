#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/9 12:16
# @Author  : Runsheng     
# @File    : flow.py

"""
Create a flow for series of NCBI SRA data
"""
# standard library import
import os
import shutil
import logging
from glob import glob

# self-import
from prep import fasta_shuffle
from bait import *


def _get_fqname(fq_dir):

    fqs=os.listdir(fq_dir)

    k_set=set()
    for fq in fqs:
        k_set.add(fq.split("_")[0])


    for fq in fqs:
        for k_single in k_set:
            if k_single in fq:
                pass


def flow_bait(i, work_dir, fq_dict, ref_file):
    """
    The flow for baiting,
    fq_dict is {lib1:[fq_1, fq_2], lib2:[fq_1]}
    todo: get
    """
    print("start")
    os.chdir(work_dir)
    ##### for each round
    xi=get_xi(i)
    print("--------"+"round"+xi+"----------------------")
    round_root_dir=os.path.join(work_dir, ("round"+xi))
    round_ref_dir=os.path.join(round_root_dir, "ref")

    round_ref_file=os.path.join(round_ref_dir, "ref.fasta")

    if os.path.exists(round_ref_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(round_ref_dir)

    os.chdir(round_root_dir)

    fasta_shuffle(ref_file=ref_file,out_file=round_ref_file)
    bwa_index_wrapper(round_ref_file)

    for k, v in fq_dict.iteritems():

        out_bam=os.path.join(round_root_dir, k+".bam")
        fq_str=" ".join(v)
        bwa_mem_wrapper(round_ref_file, fq_str, core=32, min_seed_length=18, band_width=2000, out=out_bam)
        bam_sorted=sort_index_wrapper(out_bam)
        names_set=get_name(bam_sorted)


        for fqin in v:
            fq_out=fqin.split("/")[-1].split(".")[0]+"_s.fq"
            fqout=fq_subset_main(fqin, names_set, fq_out)


def pre_fastq_dir(sra_list, work_dir=None):
    """
    step1: get the fastq file from the sra files
    :param sra_list:
    :param work_dir:
    :return:
    """
    if work_dir==None:
        work_dir=os.getcwd()
    print(work_dir)

    os.chdir(work_dir)
    fastq_raw_dir=os.path.join(work_dir, "fastq")
    if os.path.exists(fastq_raw_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(fastq_raw_dir)

    for sra_file in sra_list:
        sra2fq_wrapper(sra_file, fastq_raw_dir)

    return work_dir


def pre_ref(i, store_dir, ref_file=None):
    
    os.chdir(store_dir)

    ref_first=os.path.join(store_dir, "round0.fasta")
    
    if ref_file==None:
        print("No ref file provided, enter de novo mode, this may take a long time to finish.")
    else:
        ref_name=fasta_shuffle(ref_file=ref_file, out_file=ref_first, file_type="fasta")
        #shutil.copyfile(ref_file, ref_first) # need to note if the ref_file is in the workdir

    return ref_name


#------------------------
def get_xi(i):
    """
    used to rename the file in dir
    :param i:
    :return:
    """
    return i/10*"x"+str(i%10)


def get_i(xi):
    i=0
    for x_single in xi:
        if x_single=="x":
            i+=10
        else:
            i+=int(x_single)
    return i
#------------------------


def pre_mapping_dir(i,work_dir=None,ref_file=None):
    """
    prepare the dir and files for the i round mapping and fill
    :param i: the loop number for N_fill
    :param ref_file: the genome fasta file to be filled
    :param work_dir:
    :return: the path to work dir
    """
    if work_dir==None:
        work_dir=os.getcwd()
    print(work_dir)

    xi=get_xi(i)
    tmp_dir=os.path.join(work_dir, xi)

    tmp_ref_dir=os.path.join(tmp_dir, "ref")

    if os.path.exists(tmp_ref_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(tmp_ref_dir)

    if i==0:
        if ref_file==None:
            raise IOError("A ref file have to be provided in the first round.")
        else:
            shutil.copyfile(ref_file, (tmp_ref_dir+"/round{}.fasta".format(i)) )
    else:
        ref_last=sorted(glob((work_dir+"/scafo.fasta")))[-1]
        shutil.copy(ref_last, tmp_ref_dir)
        print(ref_last)
    return tmp_ref_dir


if __name__=="__main__":
    i=0
    work_dir="/home/zhaolab1/data/mitosra/dna/mitovar_test"
    fq_dict={"ERR1018617":["/home/zhaolab1/data/mitosra/dna/onewkdir/1094327/fastq/ERR1018617_1.fastq",
             "/home/zhaolab1/data/mitosra/dna/onewkdir/1094327/fastq/ERR1018617_2.fastq"]}
    ref_file = "/home/zhaolab1/data/mitosra/dna/ref/celcbr.fa"
    flow_bait(i, work_dir, fq_dict, ref_file)
