#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/1/20 18:57
# @Author  : Runsheng     
# @File    : submission.py

"""
used to generate the files for the submission to genbank
including
1. the .fsa with header
2. the .tbl file

the files can be converted to sqn file using tbl2asn
"""

from utils import fasta2dic, dic2fasta
from anno import *
from post import scaf_filter

def flow_submission(fastafile, rrna_ref, cds_ref, spe_name, MITFIPATH=None,outfasta=None):
    reorder_fasta=flow_re_order(fastafile, MITFIPATH, outfasta)
    prefix=reorder_fasta.split(".")[0]

    trna_file=mitfi_wrapper_trna(reorder_fasta, MITFIPATH)
    bed4_trna = trnafile_parser(trna_file)
    bed4_rrna = rrna_tbl_get(reorder_fasta, rrna_ref)
    bed4_cds = cds_tbl_get(reorder_fasta, cds_ref)

    # get all intervals separately
    tbl_out = tbl_format(bed4_rrna, bed4_cds, bed4_trna)
    # get the header line
    fastaname = spe_name.replace(" ", "_")
    headstr = ">Feature\t" + fastaname + "\n"
    # merge
    head=[]
    head.append(headstr)
    tbl_out = head + tbl_out

    # the fsa and tbl have to have the same prefix
    tbl_filename=prefix + ".tbl"
    fsa_filename=prefix + ".fsa"

    with open(tbl_filename, "w") as fw:
        fw.write("".join(tbl_out))

    sub_fsa=pre_fsa(reorder_fasta, spe_name, fsa_filename) # for output

    return (tbl_filename, fsa_filename)


def pre_fsa(fastafile, spe_name, out=None):
    """
    pre the .fsa file for sqn submission
    :param fastafile:
    :return:
    """
    fa_d=fasta2dic(fastafile)
    if out is None:
        prefix = fastafile.split(".")[0]
        out = prefix + ".fsa"
    fastaname=spe_name.replace(" ", "_")

    header=">{fastaname} [organism={spe_name}] [chromosome=mt] [moltype=genomic DNA] " \
           "[gcode=5] [Topology=Circular] [Completedness=Complete] " \
           "{spe_name} mitochondrion, complete genome.".format(
        fastaname=fastaname,spe_name=spe_name)
    with open(out, "w") as fw:
        fw.write(header)
        fw.write("\n")
        fw.write(str(fa_d.values()[0].seq))
        fw.write("\n")

    return out


if __name__=="__main__":

    # exam of the sinica
    def re_run_sinica():

        wkdir="/home/zhaolab1/data/mitosra/dna/wkdir/"
        spe="281681"
        work_dir_spe=os.path.join(wkdir, spe)
        os.chdir(work_dir_spe)
        fasta_d=scaf_filter("./round0/spades_out/scaffolds.fasta")
        fastaname=spe+".fasta"
        dic2fasta(fasta_d,fastaname)

        spe_name="Caenorhabditis plicata"
        rrna_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
        cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"

        print flow_submission(fastaname, rrna_ref, cds_ref, spe_name)

    def re_run_sinica_2nd():
        wkdir="/home/zhaolab1/data/mitosra/rnaother/sinica"
        os.chdir(wkdir)
        spe_name="Caenorhabditis sinica"
        rrna_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
        cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"
        fastaname="sinica_withdloop.fa"
        print flow_submission(fastaname, rrna_ref, cds_ref, spe_name)

    pass
