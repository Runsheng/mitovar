#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/9 12:09
# @Author  : Runsheng     
# @File    : assembler.py

from utils import myexe
import os


def spades_wrapper(fq_name_dict, outdir, core=12):

    """
    :param fq_name_dict: a dict with {"pe1-1": "ERRxxxx_1_s.fq","pe1-2": "ERRxxxx_1_s.fq", "s2":"ERRxxxx_1.fq"}
    :param outdir:
    :param core:
    :return:
    """



    spades_cmd_pe="spades.py --only-assembler -t {core} \
                    --pe1-1 fq1 \
                    --pe1-2 fq2 \
                    -o {outdir}".format(
        core=core, outdir=outdir
    )


if __name__=="__main__":
    i=0
    work_dir="/home/zhaolab1/data/mitosra/dna/mitovar_test"
    fq_dict={"ERR1018617":["/home/zhaolab1/data/mitosra/dna/onewkdir/1094327/fastq/ERR1018617_1.fastq",
             "/home/zhaolab1/data/mitosra/dna/onewkdir/1094327/fastq/ERR1018617_2.fastq"]}
    fq_s_dict={"ERR1018617":["ERR1018617_1_s.fastq,ERR1018617_2_s.fastq"]}

    ref_file = "/home/zhaolab1/data/mitosra/dna/ref/celcbr.fa"
    flow_bait(i, work_dir, fq_dict, ref_file)