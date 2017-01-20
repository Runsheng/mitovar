#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/14 20:28
# @Author  : Runsheng     
# @File    : post.py

"""
The pipelines used after a round's baiting-assembly is finished
Have to chose the reference from the scaffolds.fasta to generate new reference

"""
from utils import fasta2dic


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


def scaf_filter(filename, cutoff_length=None, cutoff_coverage=None):
    """
    get the attributes from the scaf/contig name
    :return:
    """
    fa_dict=fasta2dic(filename)
    fa_dict_f={}

    # cacl the cutoff if not given
    length_l=[]
    coverage_l=[]
    for name, seq in fa_dict.iteritems():
        name_p=scaf_name_parse(name)
        length_l.append(name_p["length"])
        coverage_l.append(name_p["coverage"])
    if cutoff_length==None:
        cutoff_length=max(length_l)/10
    if cutoff_coverage==None:
        cutoff_coverage=max(coverage_l)/20

    # do the filter
    for name, seq in fa_dict.iteritems():
        name_p=scaf_name_parse(name)
        if name_p["length"]>cutoff_length and name_p["coverage"]>cutoff_coverage:
            fa_dict_f[name]=seq

    return fa_dict_f


if __name__=="__main__":
    pass