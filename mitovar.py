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
mitovar.py anno -f mtDNA.fasta -c cel_p.fa -r cel_rrna.fa -s cel
mitovar.py assemble -f cel.fa -p 32 -s cbr

version 0.0.99
            """,
            formatter_class = argparse.RawDescriptionHelpFormatter)
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
                            help="the working dir, default is the current dir")
        parser.add_argument("-f", "--fasta", help="fasta to be annotated")
        parser.add_argument("-c", "--cds_ref",
                            help="cds sequences and/or protein sequences from one closest seed mtDNA")
        parser.add_argument("-r", "--rrna_ref", help="rrna sequences from one closest seed mtDNA")
        parser.add_argument("-s", "--spe", default="un", help="the species name, if not give, set to un")


        args = parser.parse_args(sys.argv[2:])
        #print(args.fasta)
        flow_anno(wkdir=args.wkdir, fasta=args.fasta, cds_ref=args.cds_ref,
              rrna_ref=args.rrna_ref, spe_name=args.spe)

    def assemble(self):
        parser=argparse.ArgumentParser(
            description="The command to get mtDNA aseembly from fastq NGS reads and a nearby reference, "
                        "please put a ./{spe}/fastq folder containing all fastq file inside the wkdir,"
                        "for example, you wkdir is /home and your spe is cbr, then put fastq files to"
                        "/home/cbr/fastq"
                        "the fastq files with the same prefix before _ will be treated as pair-end automatically, "
                        "such as SRR1_F.fq and SRR2_R.fq,"
                        "support at most 9 separated fastq pairs."
        )


        parser.add_argument("-d", "--wkdir", default=os.getcwd(),
                            help="the working dir, default is the current dir")
        parser.add_argument("-f", "--fasta", help="reference fasta file used to bait reads from NGS reads")
        parser.add_argument("-p", "--core", default=4, help="cores used to run mapping and assembly")
        parser.add_argument("-s", "--spe", default="new", help="the species name, if not give, set to new")


        args = parser.parse_args(sys.argv[2:])

        if "fastq" not in os.listdir(args.wkdir+"/"+args.spe):
            print("Please put all fastq files in the ./fastq forlder inside the species wkdir,"
                  "for example, you wkdir is /home and your spe is cbr, then put fastq files to"
                  "/home/cbr/fastq")
            exit(1)

        flow_chain_fq_first(args.spe,args.fasta,args.wkdir,args.core,i=0,min_seed_length=50, band_width=2000)



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

