#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/12 17:30
# @Author  : Runsheng     
# @File    : __collect.py

"""
to collect
"""
import os

def __pass():
    pass


def _get_fqname(fq_dir):

    fqs=os.listdir(fq_dir)

    k_set=set()
    for fq in fqs:
        k_set.add(fq.split("_")[0])


    for fq in fqs:
        for k_single in k_set:
            if k_single in fq:
                pass
