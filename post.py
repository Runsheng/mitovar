#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/14 20:28
# @Author  : Runsheng     
# @File    : post.py

"""
The pipelines used after a round's baiting-assembly is finished
Have to chose the reference from the scaffolds.fasta to generate new reference

"""
from utils import fasta2dic, myexe
from flow import get_fq_dict
import os
from bait import *

def scaf_name_parse(name):
    """
    ">NODE_1_length_13188_cov_721.856"
    :param name:
    :return:
    """
    name_l=name.split("_")
    _, node_num, _, length, _, cov=name_l
    coverage=cov.split(".")[0]
    return {"node":node_num, "length": int(length), "coverage":int(coverage)}


def scaf_filter(filename, cutoff_length=None, cutoff_coverage=None, len_cutoff=13000):
    """
    get the attributes from the scaf/contig name
    :return:
    """
    fa_dict=fasta2dic(filename)
    fa_dict_f={}

    # cacl the cutoff if not given
    length_l=[]
    coverage_l=[]
    total_len=0
    for name, seq in fa_dict.iteritems():
        name_p=scaf_name_parse(name)
        length_l.append(name_p["length"])
        coverage_l.append(name_p["coverage"])
    if cutoff_length is None:
        cutoff_length=max(length_l)/10
    if cutoff_coverage is None:
        cutoff_coverage=max(coverage_l)/20

    # do the filter
    for name, seq in fa_dict.iteritems():
        name_p=scaf_name_parse(name)
        if name_p["length"]>cutoff_length and name_p["coverage"]>cutoff_coverage\
                and total_len<len_cutoff:
            fa_dict_f[name]=seq
            total_len+=len(seq)

    return fa_dict_f


def flow_mapping(work_dir, fq_dict, ref_file, core=16):
    """
    The flow for baiting,
    fq_dict is {lib1:[fq_1, fq_2], lib2:[fq_1]}
    todo: get
    """
    print("start")
    os.chdir(work_dir)
    ##### for each round
    print("--------"+"mapping and merge samfile"+"----------------------")

    bwa_index_wrapper(ref_file)

    bam_list=[]
    for k, v in fq_dict.items():
        fq_str = " ".join(v)
        out_bam = os.path.join(work_dir, k + ".bam")
        bwa_mem_wrapper(ref_file, fq_str,
                              core=core, min_seed_length=20,
                              band_width=2000, out=out_bam)
        bam_sorted=sort_index_wrapper(out_bam)
        bam_list.append(bam_sorted)

    return bam_list


def wrapper_samtools_merge(ref_file, bam_list, out=None):
    if out is None:
        out=ref_file.split(".")[0]+".bam"

    merge_cmd="samtools merge {out_bam} {in_bam}".format(
        out_bam=out, in_bam=" ".join(bam_list))
    print(merge_cmd)
    myexe(merge_cmd)

    return sort_index_wrapper(out)


def post_mapping(work_dir, fq_dir, ref_file, core=32 , out=None):
    """
    bwa mapping and merge all sam files for a spe
    :return:
    """
    os.chdir(work_dir)
    fq_dict=get_fq_dict(fq_dir)
    print "The fq contains:", fq_dict

    bam_list=flow_mapping(work_dir, fq_dict, ref_file, core)

    out=wrapper_samtools_merge(ref_file, bam_list, out)
    return out



if __name__=="__main__":
    dir_list = ['281687', '1561998', '1729975', '1094327', '1094328',
                '1094335', '1094320', '1094321', '1094326', '135651', '860376']
    for taxid in dir_list:
        ref_file="/home/zhaolab1/data/mitosra/dna/ref_post/{taxid}_s_ordered.fsa".format(taxid=taxid)
        work_dir_spe="/home/zhaolab1/data/mitosra/dna/wkdir/{taxid}".format(taxid=taxid)
        print post_mapping(work_dir_spe, ref_file, core=40)