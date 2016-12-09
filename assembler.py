#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/9 12:09
# @Author  : Runsheng     
# @File    : assembler.py

from utils import myexe
import os


def spades_wrapper(fq_name_dict, outdir, core=12):

    spades_cmd_pe="spades.py --only-assembler -t {core} \
                    --pe1-1 fq1 \
                    --pe1-2 fq2"