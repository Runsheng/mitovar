#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/9 13:08
# @Author  : Runsheng     
# @File    : anno.py
"""
Try to annotate the assembly of mito scafs, including the tRNA, rRNA and coding genes
Find the tRNA-phe,  set it as the start point of the new reference genome
"""
from utils import fasta2dic, chr_select, myexe, dic2fasta
from Bio import SearchIO, SeqIO, pairwise2
import os

def exonerate_wrapper(query, target, outfile=None, geneticcode=5):
    """
    --geneticcode 5

    return is a outfile name in relative path
    todo: using stringIO to hinder the file IO
    """
    if outfile is None:
        outfile=query.split("/")[-1].split(".")[0]+".exonerate"

    exonerate_cmd="exonerate {query} {target} \
                   --geneticcode {geneticcode} \
                   > {outfile}".format(
                        query=query, target=target,
                        geneticcode=geneticcode,
                        outfile=outfile)
    myexe(exonerate_cmd)

    return outfile


def exonerate_parser_write(query, exonerate_file, prefix=None):
    """
    parser the exonerate result, and return the protein and cds file
    modification: add the position to the name of the cds and pro, use space to add interval
    :param query:
    :param exonerate_file:
    :param prefix:
    :return:
    """
    ref_dict=fasta2dic(query)

    if prefix is None:
        prefix=exonerate_file.split(".")[0]

    p_outname=(prefix+"_exonerate_p.fa")
    cds_outname=(prefix + "_exonerate_cds.fa")

    fw_p=open(p_outname, "w")
    fw_cds=open(cds_outname, "w")

    texts=SearchIO.parse(exonerate_file, format="exonerate-text")
    for record in texts:
        for hsp in record:
            for s in hsp:
                #print(s.fragment.hit_id)
                name_str=">"+s.fragment.hit_id
                name_cds, cds_str=chr_select(ref_dict, s.fragment.query_id, s.fragment.query_start,
                                 s.fragment.query_end)
                p_str=str(s.fragment.query.seq)

                #print(name_cds)

                fw_p.write(name_str+"  "+name_cds+"\n"+p_str+"\n")
                fw_cds.write(name_str+"  "+name_cds+"\n"+cds_str+"\n")

    return cds_outname, p_outname


def exon_corrector(cds_outname, p_outname, prefix=None):
    """
    try to merge the cds and the protein sequence together,
    if, for example, the ND3 is divided into two part,
    The easiest way is to just add the two segment together
    :param cds_outname:
    :param p_outname:
    :return:
    todo: unfinished, only chose the longest sequence
    """
    cds_d={}
    p_d={}
    if prefix is None:
        prefix=cds_outname.split("_")[0]

    for cds, p in zip(SeqIO.parse(cds_outname, "fasta"), SeqIO.parse(p_outname, "fasta")):
        name, anno=cds.description.split("  ")
        __=anno.split(":") # unused
        key="#".join([prefix, name])
        if key in cds_d.keys():
            seq_old=cds_d[key]
            if len(cds)>len(seq_old):
                cds_d[key]=cds
                p_d[key]=p
            else:
                pass
        else:
            cds_d[key]=cds
            p_d[key]=p

    dic2fasta(cds_d, "_".join([prefix, "cds", "corr.fasta"]))
    dic2fasta(p_d, "_".join([prefix, "p", "corr.fasta"]))


def genbank_parser(gbfile, prefix=None):
    """
    from gb file, get the cds and protein sequence
    and rename them as #
    :param gbfile: genbank file
    :return:
    """
    cds_d={}
    p_d={}

    if prefix is None:
        prefix = gbfile.split(".")[0]

    for record in SeqIO.parse(gbfile, format="genbank"):
        for i in record.features:
            if i.type == "CDS":
                name = "".join(i.qualifiers["gene"])
                p = "".join(i.qualifiers['translation'])
                cds = str(i.extract(record.seq))

                key = "#".join([prefix, name])
                if key in cds_d.keys():
                    seq_old = cds_d[key]
                    if len(cds) > len(seq_old):
                        cds_d[key] = cds
                        p_d[key] = p
                    else:
                        pass
                else:
                    cds_d[key] = cds
                    p_d[key] = p

    dic2fasta(cds_d, "_".join([prefix, "cds", "corr.fasta"]))
    dic2fasta(p_d, "_".join([prefix, "p", "corr.fasta"]))


def flow_exon(query, target, outfile=None, geneticcode=5, prefix=None):
    """
    get the CDS region and translated protein after using,
    the default geneticcode is 5, the invertebrate mitochondrial code
    :param query:
    :param target:
    :param outfile:
    :param geneticcode:
    :param prefix:
    :return:
    """
    outfile=exonerate_wrapper(query, target, outfile, geneticcode)
    cds_out, p_out=exonerate_parser_write(query, outfile, prefix)
    exon_corrector(cds_out, p_out)


def get_cogfile(fastafile, wkdir=None, out="m20.txt"):
    """
    cat all the cds/or protein together, and get the cog file used for ete alignment
    :param fastafile: the comibined fasta file
    :param out: the name indicaiting the orthologs, like: "ppac#ND4\tcele#ND4\nppac#ND5\tcele#ND5\n"
    :return:
    """
    if wkdir is None:
        wkdir=os.getcwd()

    os.chdir(wkdir)
    fa_d=fasta2dic(fastafile)
    fw=open(out, "w")
    name_d={}
    for k in fa_d.keys():
        suffix=k.split("#")[1]
        try:
            name_d[suffix].append(k)
        except KeyError:
            name_d[suffix]=[]
            name_d[suffix].append(k)
    for k, v in name_d.iteritems():
        fw.write("\t".join(v))
        fw.write("\n")

    return out


if __name__=="__main__":

    def test_flowexon():
        from glob import glob
        import os
        wkdir="/home/zhaolab1/data/mitosra/dna/anno/exon"
        os.chdir(wkdir)
        target="./ref/celmt_p.fasta"
        query_list=glob("*_s.fasta")
        print(query_list)
        for query in query_list:
            print query
            flow_exon(query, target)
        myexe("""
        cat *_p_corr.fasta> m20_p.fasta
        cat *_cds_corr.fasta >m20_cds.fasta
        """)

    def test_corr():
        wkdir="/home/zhaolab1/data/mitosra/dna/anno/exon/round0"
        os.chdir(wkdir)
        exon_corrector("281687_s_cds.fa","281687_s_p.fa")

    #test_flowexon()
    def test_cog():
        get_cogfile("m20_p.fasta")

    def test_genbank_parser():
        wkdir="/home/zhaolab1/data/mitosra/dna/anno/exon"
        os.chdir(wkdir)


    test_flowexon()
