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
from utils import myglob, fasta2dic, dic2fasta

# self-import
from prep import fasta_shuffle
from bait import *
from assembler import spades_wrapper
from anno import exonerate_parser_write, exonerate_wrapper
from post import scaf_filter

def flow_bait(i, work_dir, fq_dict, ref_file, core=16):
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

    fqout_dict={} # a dict with {prefix1:[read1, read2, read3], prefix2:[read1]},
    for k, v in fq_dict.iteritems():

        out_bam=os.path.join(round_root_dir, k+".bam")
        fq_str=" ".join(v)
        bwa_mem_wrapper(round_ref_file, fq_str, core=core, min_seed_length=20, band_width=2000, out=out_bam)
        bam_sorted=sort_index_wrapper(out_bam)
        names_set=get_name(bam_sorted)

        fq_pair=[]
        for fqin in v:
            fq_out_name=fqin.split("/")[-1].split(".")[0]+"_s.fq" # fqin could be abs path like /home/user/read1.fq
            fq_subset_main(fqin, names_set, fq_out_name)
            fq_pair.append(fq_out_name)
        fqout_dict[k]=fq_pair
    return fqout_dict


def flow_pre_fastq(sra_list, sra_dir, work_dir=None):
    """
    step0: get the fastq file from the sra files
    :param sra_list: a list of sra file for one species, [ERRxxxx, SRRxxxx...]
    :param sra_dir: the dir storing the sra file, if blank, then try to using the online sra database,
                    Warning, using the online sra file can greatly increase the runtime!
    :param work_dir:
    :return: fq_dict: a dict with {"ERRxxxx": ["ERRxxxx_1.fq","ERRxxxx_2.fq"], "ERRxxx2":["ERRxxx2_1.fq"]}
    """
    if work_dir==None:
        work_dir=os.getcwd()
    print(work_dir)
    sra_file_l = [os.path.join(sra_dir, x) for x in sra_list]

    os.chdir(work_dir)
    fastq_raw_dir=os.path.join(work_dir, "fastq")
    if os.path.exists(fastq_raw_dir):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(fastq_raw_dir)
    for sra_file in sra_file_l:
        sra2fq_wrapper(sra_file, fastq_raw_dir)
    fq_dict=get_fq_dict(fastq_raw_dir)

    return fq_dict


def flow_chain_sra_scaf(spe, sra_list, ref_file, work_dir_root, sra_dir="", core=16, i=0):
    """
    The first loop of the protocol, from the sra files
    :param spe:
    :param sra_list:
    :param ref_file:
    :param work_dir_root:
    :param sra_dir:
    :param core:
    :param i:
    :return:
    """
    work_dir_spe=os.path.join(work_dir_root, spe)
    if os.path.exists(work_dir_spe):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(work_dir_spe)

    fq_dict=flow_pre_fastq(sra_list, sra_dir, work_dir_spe)

    fq_out_dict=flow_bait(i, work_dir_spe, fq_dict, ref_file)
    scaf_fasta, scaf_fastg=spades_wrapper(fq_name_dict=fq_out_dict, core=core, outdir="spades_out")
    return scaf_fasta, scaf_fastg


def pre_ref_first(store_dir, ref_file=None):
    """
    only used for round0
    :param store_dir:
    :param ref_file:
    :return:
    """
    os.chdir(store_dir)

    ref_first=os.path.join(store_dir, "round0.fasta")
    
    if ref_file==None:
        print("No ref file provided, enter de novo mode, this may take a long time to finish.")
    else:
        ref_name=fasta_shuffle(ref_file=ref_file, out_file=ref_first, file_type="fasta")
        #shutil.copyfile(ref_file, ref_first) # need to note if the ref_file is in the workdir

    return ref_name


def pre_ref_next(i, work_dir_spe):
    """
    join the reference from the last time and this time's scaffold together
    expect a dir str like:
    work_dir
        -round(i-1)
            -ref
                -roundi-2.fasta
            -spades_out
                -scaffold.fasta
    :param i:
    :param work_dir_spe:
    :return:
    """
    i_xi=get_xi(i)
    work_dir_round=os.path.join(work_dir_spe, ("round" + i_xi))

    os.chdir(work_dir_round)
    print(work_dir_round)

    ref_old=myglob(work_dir_round, "ref.fasta")[0]
    ref_new=myglob(work_dir_round, "scaffolds.fasta")[0]
    print ref_new

    ref_new_d=scaf_filter(ref_new)
    ref_old_d=fasta2dic(ref_old)

    len_scaf=0
    for v in ref_new_d.values():
        len_scaf+=len(v.seq)
    if len_scaf<12000.0: # only include the old ref when the generated scaf is short
        print("Using mixed old-new reference!")
        ref_new_d.update(ref_old_d)

    out=os.path.join(work_dir_round, "ref_next.fasta")
    dic2fasta(ref_new_d, out)

    return out


def flow_chain_general(spe, work_dir_root,ref_file=None, core=32, i_start=1, i_end=9):
    """
    The other loop of the protocol,
    the str of the work_dir_root is expected as:
        - fastq
        - round0
        - roundi

    loop of the protocol, from the sra files

    retrun: the last round of the scaf, and the fastg file for the assembly
    """
    work_dir_spe=os.path.join(work_dir_root, spe)
    os.chdir(work_dir_spe)

    # get the round0 ref
    if ref_file is None:
        ref_next=pre_ref_next(i_start-1, work_dir_spe)
    else:
        ref_next=ref_file
    work_dir_fq=os.path.join(work_dir_spe, "fastq")

    for i in range(i_start, i_end):
        fq_dict=get_fq_dict(work_dir_fq)
        fq_out_dict=flow_bait(i, work_dir_spe, fq_dict, ref_next, core=core)
        scaf_fasta, scaf_fastg=spades_wrapper(fq_name_dict=fq_out_dict, core=core, outdir="spades_out")
        ref_next=pre_ref_next(i, work_dir_spe)

    return scaf_fasta, scaf_fastg


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
    def test_1_round0():
        i=0
        work_dir="/home/zhaolab1/data/mitosra/dna/mitovar_test"
        fq_dict={"ERR1018617":["/home/zhaolab1/data/mitosra/dna/onewkdir/1094327/fastq/ERR1018617_1.fastq",
                 "/home/zhaolab1/data/mitosra/dna/onewkdir/1094327/fastq/ERR1018617_2.fastq"]}
        ref_file = "/home/zhaolab1/data/mitosra/dna/ref/celcbr.fa"
        flow_bait(i, work_dir, fq_dict, ref_file)
    # test round1-9
    spe="497829"
    work_dir_root="/home/zhaolab1/data/mitosra/dna/wkdir/"
    flow_chain_general(spe=spe, work_dir_root=work_dir_root, core=34, i_start=1, i_end=4)

