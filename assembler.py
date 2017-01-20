#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/9 12:09
# @Author  : Runsheng     
# @File    : assembler.py

from utils import myexe
import os


def spades_wrapper(fq_name_dict, outdir="spades_out", core=12):

    """
    :param fq_name_dict: a dict with {"pe1-1": "ERRxxxx_1_s.fq","pe1-2": "ERRxxxx_1_s.fq", "s2":"ERRxxxx_1.fq"}
    :param outdir:
    :param core:
    :return: the file position of the scaf fasta and the fastg file
    """
    fq_str_list=[]
    # to generate the readpool string used for spades
    lib_num=1
    for k, v in fq_name_dict.iteritems():
        read_type="pe" if len(v)>=2 else "s"
        for fq_one_name in v:
            if read_type=="s":
                str_one="--{read_type}{lib_num} {fq_one_name}".format(
                        read_type=read_type, lib_num=lib_num, fq_one_name=fq_one_name)
            elif read_type=="pe":
                if "_1" in fq_one_name:
                    fq_pos_ind="-1"
                elif "_2" in fq_one_name:
                    fq_pos_ind="-2"
                elif "_3" in fq_one_name:
                    fq_pos_ind="-s"

                str_one="--{read_type}{lib_num}{fq_pos_ind} {fq_one_name}".format(
                        read_type=read_type, lib_num=lib_num, fq_pos_ind=fq_pos_ind, fq_one_name=fq_one_name)
            fq_str_list.append(str_one)

        lib_num+=1
        fq_str=" ".join(fq_str_list)

    spades_cmd="spades.py --careful -t {core} {readpool} -o {outdir}".format(
                    core=core, readpool=fq_str, outdir=outdir)
    print(spades_cmd)
    print(myexe(spades_cmd))

    scaf_fasta=os.path.join(outdir, "scaffolds.fasta")
    scaf_fastg=os.path.join(outdir, "assembly_graph.fastg")
    return scaf_fasta, scaf_fastg


if __name__=="__main__":
    pass
    work_dir="/home/zhaolab1/data/mitosra/dna/mitovar_test/round0"
    os.chdir(work_dir)
    fq_out_dict={"ERR1018617":["/home/zhaolab1/data/mitosra/dna/mitovar_test/round0/ERR1018617_1_s.fq",
             "/home/zhaolab1/data/mitosra/dna/mitovar_test/round0/ERR1018617_2_s.fq"]}
    work_dir="/home/zhaolab1/data/mitosra/dna/wkdir/1094320"

    spades_wrapper(fq_name_dict=fq_out_dict, core=32, outdir="spades_out")
