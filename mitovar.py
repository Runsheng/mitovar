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
from flow import flow_chain_fq_first


def flow_anno(wkdir, fasta, cds_ref, rrna_ref, spe_name):
    os.chdir(wkdir)
    out=flow_re_order(fastafile=fasta)
    flow_tbl(out, cds_ref, rrna_ref, spe_name)


class CMD(object):

    def __init__(self):
        parser=argparse.ArgumentParser(
            description="Mitovar cmd lines",
            usage=""" mitovar <command> [<args>]

-------
The command contains:
anno: annotate a mtDNA fasta sequences and generate a tbl file for genbank submission
assemble: assemble the mtDNA fasta file from a NGS fastq file and a nearby reference file
------
A example for running annotation command:
mitovar anno -f mtDNA.fasta
            """,
            formatter_class = argparse.RawTextHelpFormatter)
        parser.add_argument("command", help="Subcommand to run")
        args=parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)
        getattr(self, args.command)()


    def anno(self):
        parser=argparse.ArgumentParser(
            description="The command to annotate the mtDNA fasta file"
        )


        parser.add_argument("-d", "--wkdir", default=os.getcwd(),
                            help="the working dir, if not give, set to current dir")
        parser.add_argument("-f", "--fasta", help="fasta to be annotated")
        parser.add_argument("-c", "--cds_ref",
                            help="cds sequences and/or protein sequences from one closest seed mtDNA")
        parser.add_argument("-r", "--rrna_ref", help="rrna sequences from one closest seed mtDNA")
        parser.add_argument("-p", "--spe", default="un", help="the species name, if not give, set to un")



        args = parser.parse_args(sys.argv[2:])
        #print(args.fasta)
        flow_anno(wkdir=args.wkdir, fasta=args.fasta, cds_ref=args.cds_ref,
              rrna_ref=args.rrna_ref, spe_name=args.spe)


    def assemble(self):
        pass



def __main():
    print("---The cmd wrapper for mitovar package---")
    parent_parser=argparse.ArgumentParser(prog="mitovar")
    parent_parser.add_argument("--version", action="store_true", help="mitovar version 0.0.99")

    parser = argparse.ArgumentParser(add_help=False)
    subparsers= parent_parser.add_subparsers(help=False)

    parser_anno= subparsers.add_parser("anno",  parents=[parent_parser])

    parser_anno.add_argument("-d", "--wkdir", default=os.getcwd(),
                             help="the working dir, if not give, set to current dir")
    parser_anno.add_argument("-f","--fasta", help="fasta to be annotated")
    parser_anno.add_argument("-c","--cds_ref",
                             help="cds sequences and/or protein sequences from one closest seed mtDNA")
    parser_anno.add_argument("-r", "--rrna_ref", help="rrna sequences from one closest seed mtDNA")
    parser_anno.add_argument("-p", "--spe", default="un", help="the species name, if not give, set to un")

    args=parser.parse_args()
    #args= parser_anno.parse_args()
    flow_anno(wkdir=args.wkdir, fasta=args.fasta, cds_ref=args.cds_ref,
              rrna_ref=args.rrna_ref, spe_name=args.spe)



if __name__=="__main__":
    CMD()

    ##### example for anno
   # cd home/zhaolab1/data/mitosra/dna/wkdir/1094320/round0

   # /home/zhaolab1/myapp/mitovar/mitovar.py anno \
   #     --wkdir="/home/zhaolab1/data/mitosra/dna/wkdir/1094320/round0" \
   #     --r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta" \
   #     --cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta" \
   #     --fasta="1094320.fasta" \
   #     --spe_name="sp1"
    #####

