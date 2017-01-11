#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/8 17:01
# @Author  : Runsheng     
# @File    : prep.py

import os
import shutil
from Bio import SeqIO
from utils import myexe
from glob import glob


def check_depend():
    """
    todo: check the required software needed for mitovar
    :return:
    """
    pass


def fasta_shuffle(ref_file, out_file, file_type="fasta"):

    ref_dict=SeqIO.to_dict(SeqIO.parse(ref_file, file_type))

    fw=open(out_file, "w")

    k_used_set=set()
    for k in ref_dict.keys():
        if "shuffled" in k:
            k_used_set.add(k.replace("_shuffled", ""))

    for k,v in ref_dict.iteritems():

        if "shuffled" in k or k in k_used_set:
            fw.write(">" + k + "\n")
            fw.write(str(v.seq))
            fw.write("\n")
        else:
            v_str=str(v.seq)
            v_new=v_str[len(v)/2:]+v_str[:len(v)/2]
            k_new=str(k)+"_shuffled"

            fw.write(">"+k+"\n")
            fw.write(v_str)
            fw.write("\n")

            fw.write(">"+k_new+"\n")
            fw.write(v_new)
            fw.write("\n")

    fw.close()

    return out_file


if __name__=="__main__":
    pass