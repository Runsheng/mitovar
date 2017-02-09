#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/12 14:50
# @Author  : Runsheng     
# @File    : flow_test.py

from flow import flow_bait
from assembler import spades_wrapper
import os
from flow import flow_pre_fastq
from flow import flow_chain_sra_scaf, flow_chain_fq_first
from bait import *
from collections import OrderedDict
from anno import flow_exon


def test1():
    # just test the reads from previous combines, treat them as single
    # and do one round of spades assembler
    work_dir_root="/home/zhaolab1/data/mitosra/dna/wkdir"
    dir_list=['281687', '1561998', '1729975', '1094327', '1094328', '31234',
             '1094335', '497829', '1094320', '1094321', '1094326', '135651', '860376']
    #dir_list=['860376']
    ref_file = "/home/zhaolab1/data/mitosra/dna/ref/celcbr.fa"
    for spe in dir_list:
        work_dir_spe_root=os.path.join(work_dir_root, spe)
        os.chdir(work_dir_spe_root)
        fq_dict={"all":["/home/zhaolab1/data/mitosra/dna/wkdir/{}/all.fastq".format(spe)]}
        fq_out_dict=flow_bait(0, work_dir=work_dir_spe_root, fq_dict=fq_dict, ref_file=ref_file)

        work_dir_round="/home/zhaolab1/data/mitosra/dna/wkdir/{}/round0".format(spe)
        os.chdir(work_dir_round)
        spades_wrapper(fq_name_dict=fq_out_dict, core=32, outdir="spades_out_onlyas")


def run_single():
    # just test the reads from previous combines, treat them as single
    # and do one round of spades assembler
    work_dir_root="/home/zhaolab1/data/mitosra/dna/wkdir"
    #dir_list=['281687', '1561998', '1729975', '1094327', '1094328', '31234',
    #         '1094335', '497829', '1094320', '1094321', '1094326', '135651', '860376']
    dir_list=['281687']
    ref_file = "/home/zhaolab1/data/mitosra/dna/ref/merge.fasta"
    for spe in dir_list:
        flow_chain_fq_first(spe, ref_file, work_dir_root, core=32)


def test_pre_fastq():
    """
    using 31243, remani genome to get the fastq
    :return:
    """
    sralist=['SRR1056292', 'SRR1056294', 'SRR1056295', 'SRR1056293', 'SRR275642']
    work_dir_root="/home/zhaolab1/data/mitosra/dna/wkdir/31234"
    flow_pre_fastq(sralist, work_dir_root)


def test_as_one():
    #spe="860376"
    spe_l=["281681","281687","497829","1094320","1094327","1094331"]
    for spe in spe_l:
        work_dir_root="/home/zhaolab1/data/mitosra/dna/wkdir/"
        ref_file = "/home/zhaolab1/data/mitosra/dna/ref/merge.fasta"
        print flow_chain_fq_first(spe, ref_file, work_dir_root, sra_dir="", core=32, i=0)


def test_half():
    dict_half=OrderedDict([
                 ('1729975', ['ERR1055244.sra']),
                 ('1094331',
                  ['ERR1059219.sra',
                   'ERR1059222.sra',
                   'ERR1059220.sra',
                   'ERR1059221.sra']),
                 ('281681',
                  ['ERR1059191.sra', 'ERR1059192.sra', 'ERR1059193.sra'])])
    sra_dir = "/home/zhaolab1/data/mitosra/sra/mitodna"
    work_dir_root="/home/zhaolab1/data/mitosra/dna/wkdir"
    ref_file = "/home/zhaolab1/data/mitosra/dna/ref/celcbr.fa"

    for spe, sra_list in dict_half.iteritems():
        print(flow_chain_sra_scaf(spe, sra_list, ref_file, work_dir_root, sra_dir, core=16))


def test_round1():
    pass


def flow_round0_all():
    from glob import glob
    import os
    wkdir="/home/zhaolab1/data/mitosra/dna/anno/exon"
    os.chdir(wkdir)
    target="./ref/celmt_p.fasta"
    query_list=glob("*.fasta")

    for query in query_list:
        print query
        flow_exon(query, target)

from post import post_mapping

def test_post_mapping():
    dir_list =['860376', '1729975', '1094326', '1502938',
               '281681', '1094328', '1094335', '1094320',
               '1094321', '1094331', '1094327']
    for taxid in dir_list:
        ref_file="/home/zhaolab1/data/mitosra/rnaother/ref/{taxid}.fasta".format(taxid=taxid)
        work_dir_spe="/home/zhaolab1/data/mitosra/rnaother/{taxid}".format(taxid=taxid)
        fq_dir=work_dir_spe+"/fastq"
        print post_mapping(work_dir_spe, fq_dir, ref_file, core=40)

def run_cbr_mapping():
    fqdir="/home/zhaolab1/data/mitosra/rna/cbr_cni/read/briggsae_male/fastq/"
    ref_file="/home/zhaolab1/data/mitosra/rna/cbr_cni/ref/cbrmt.fa"
    wkdir="/home/zhaolab1/data/mitosra/rna/cbr_cni"
    print post_mapping(work_dir=wkdir, fq_dir=fqdir, ref_file=ref_file, core=40)

def run_cni_mapping():
    fqdir="/home/zhaolab1/data/mitosra/rna/cbr_cni/read/nigoni_male/fastq/"
    ref_file="/home/zhaolab1/data/mitosra/rna/cbr_cni/ref/cnimt.fa"
    wkdir="/home/zhaolab1/data/mitosra/rna/cbr_cni"
    print post_mapping(work_dir=wkdir, fq_dir=fqdir, ref_file=ref_file, core=40)


if __name__=="__main__":
    run_cni_mapping()
    # note 86 should have another round to finish