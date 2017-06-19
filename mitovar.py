#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/5/17 12:58
# @Author  : Runsheng
# @File    : mitovar.py


"""
The runfile containing the executable scripts to
- assemble the mtDNA sequences from fastq read files or from the sra files
    -output:
        fasta file for the assembly
        fastg file to show the final graph
- annotate a near-finished fasta file
        adjust the start point to tRNA-p
        vcf file for each libraries
        tbl file for the annotation
    required field: wkdir, fasta, rrna_ref, cds_ref

- resume last run
"""

import argparse
import sys
import os

from anno import flow_re_order, flow_tbl

def flow_anno(wkdir, fasta, cds_ref, rrna_ref, spe_name):
    os.chdir(wkdir)
    out=flow_re_order(fastafile=fasta)
    flow_tbl(out, cds_ref, rrna_ref, spe_name)

if __name__=="__main__":
    import os

    parser=argparse.ArgumentParser(prog="mitovar")
    parser.add_argument("--version", action="store_true", help="mitovar version 0.0.99")
    subparsers= parser.add_subparsers(help="sub-command as anno, var and assemble")

    parser_anno= subparsers.add_parser("anno", help="annotation the mtDNA assembly")
    parser_anno.add_argument("wkdir", default=os.getcwd(), help="the working dir, if not give, set to current dir")
    parser_anno.add_argument("fasta", help="fasta to be annotated")
    parser_anno.add_argument("cds_ref", help="cds sequences from one closest seed mtDNA")
    parser_anno.add_argument("rrna_ref", help="rrna sequences from one closest seed mtDNA")
    parser_anno.add_argument("spe", default="un", help="the species name, if not give, set to un")

    flow_anno(wkdir)


    wkdir="/home/zhaolab1/data/mitosra/dna/wkdir/1094320/round0"

    os.chdir()
    r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
    cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"
    out = flow_re_order("1094320.fasta")
    flow_tbl(out, cds_ref, r_ref, spe_name="Caenorhabditis sp.1")