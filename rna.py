#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/2/22 12:44
# @Author  : Runsheng     
# @File    : rna.py

"""
Collections of functions and flows used in the RNA-seq data
The functions used to assemble the CDS sequences using RNA-seq data

"""
import os

from assembler import spades_wrapper
from flow import get_fq_dict, flow_bait


def flow_chain_rna(spe, ref_file, work_dir_root, sra_dir="", core=16, i=0):

    work_dir_spe=os.path.join(work_dir_root, spe)
    if os.path.exists(work_dir_spe):
        print("Already have the dirs.")
        pass
    else:
        os.makedirs(work_dir_spe)
    workdir_fastq=os.path.join(work_dir_spe, "fastq")

    fq_dict=get_fq_dict(workdir_fastq)
    fq_out_dict=flow_bait(i, work_dir_spe, fq_dict, ref_file, core=core)
    scaf_fasta, scaf_fastg=spades_wrapper(fq_name_dict=fq_out_dict, core=core, outdir="spades_out", rna_model=True)

    return scaf_fasta, scaf_fastg



if __name__=="__main__":
    pass