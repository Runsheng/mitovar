#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/9 13:08
# @Author  : Runsheng     
# @File    : anno.py
"""
Try to annotate the assembly of mito scafs, including the tRNA, rRNA and coding genes
Find the tRNA-phe,  set it as the start point of the new reference genome.
Get the annotation for rRNA and tRNA, store as .tbl file, used for tbl2asn submission
Get the annotation for CDS and proteins, store as .tbl file, used for tbl2asn submission.
Write the CDS and protein sequence out as fasta file, together with the orthology table, used for tree construction.
"""
from utils import fasta2dic, chr_select, myexe, dic2fasta

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO, SeqIO, SeqUtils


import os
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO


##################CDS part
def exonerate_wrapper(query, target, outfile=False, geneticcode=5, score=100, bestn=None):
    """
    --geneticcode 5

    return is a outfile name in relative path
    todo: using stringIO to hinder the file IO
    """
    if bestn is None:
        bestn=len(fasta2dic(target)) # default, output one region for one query

    exonerate_cmd="exonerate {query} {target} \
                   --geneticcode {geneticcode} \
                   --score {score} \
                   --bestn {bestn} \
                   ".format(
                        query=query, target=target,
                        geneticcode=geneticcode,
                        score=score,
                        bestn=bestn,
                        )
    out=myexe(exonerate_cmd)

    ## trigger to write the outfile to disk
    if outfile:
        outname=query.split("/")[-1].split(".")[0]+".exonerate"
        with open(outname, "w") as fw:
            fw.write(outname)

    return out


def exonerate_parser(exonerate_file):
    """
    parser the exonerate result, and return the position of the feather in 4-col bed format
    4 col bed4: [chro, start,end, name], example ["seq1", 1, 55, "trnP"]
    :param query:
    :param exonerate_file:
    :param prefix:
    :return: list of bed4
    """
    #fw=open(tbl_outname, "w") # change IO to list store
    bed4=[]

    texts=SearchIO.parse(StringIO(exonerate_file), format="exonerate-text")
    try:
        for record in texts:
            for hsp in record:
                for s in hsp:
                    # the biopython.SearchIO interval is 0 based [start, end), so start+1, end+0 to get 1 based coords
                    table_4=[s.fragment.query_id, s.fragment.query_start+1, s.fragment.query_end,s.fragment.hit_id]
                    bed4.append(table_4)
    except ValueError as e:
        pass
    bed4.sort()
    return bed4


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
        prefix=query.split(".")[0]

    p_outname=(prefix+"_exonerate_p.fa")
    cds_outname=(prefix + "_exonerate_cds.fa")

    fw_p=open(p_outname, "w")
    fw_cds=open(cds_outname, "w")

    texts=SearchIO.parse(StringIO(exonerate_file), format="exonerate-text")
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


def _exon_corrector(cds_outname, p_outname, prefix=None):
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
    #_exon_corrector(cds_out, p_out) # in default model, no need to select again

####rrna

def rrna_tbl_get(query, target, outfile=False, geneticcode=5):
    """
    rrna bed 4
    """
    out=exonerate_wrapper(query, target, outfile, geneticcode)
    tbl_rrna=exonerate_parser(out)
    return tbl_rrna


def cds_tbl_get(query, target, outfile=False, geneticcode=5):
    out= exonerate_wrapper(query, target, outfile, geneticcode)
    tbl_cds = exonerate_parser(out)
    return tbl_cds


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
#################### end of CDS part

#### end of rrna


#################### trna part
def mitfi_wrapper_trna(fastafile, MITFIPATH=None, prefix=None):
    """
    mitfi.jar in in $MITFIPATH=./bins
    :return:teh filename of mitfi run
    """
    if MITFIPATH is None:
        path = os.path.dirname(__file__)
        MITFIPATH=os.path.join(path, "bins", "mitfi")
        #print MITFIPATH
    jarfile=os.path.join(MITFIPATH, "mitfi.jar")

    mitfi_cmd="java -jar {jarfile} {fastafile}".format(
        jarfile=jarfile, fastafile=fastafile)
    trna_out=myexe(mitfi_cmd)

    print(trna_out)
    if prefix is None:
        prefix=".".join(fastafile.split("/")[-1].split(".")[0:-1])
    with open(prefix+"_trna.txt", "w") as fw:
        fw.write(trna_out)

    return prefix+"_trna.txt"


def trnafile_parser(trnafile):
    bed4_l=[]
    with open(trnafile, "r") as fr:
        for line in fr.readlines():
            if line.startswith("#"):
                pass
            else:
                terms = line.strip().split("\t")
                header, start, stop, score, evalue, AC, AA, model, strand = terms
                seq1=AA[0]
                bed4_l.append([header, int(start), int(stop), seq1])
    return bed4_l


def tbl_format(bed4_rrna, bed4_cds, bed4_trna):
    """
    tbl format :
    ---
    >refname # once
    ---
    for each term: 2line anntation
    start\tend\ttype\n\t\t\tkey\tvalue\n
    ---
    trna and rrna shows once,
    but cds show as gene and cds

    :param bed4_rrna:
    :param bed4_cds:
    :param bed4_trna:
    :return:
    """
    #sanity check
    if bed4_rrna[0][0]==bed4_cds[0][0]==bed4_trna[0][0]:
        ref=bed4_rrna[0][0]
    else:
        return "Error, annotations not from the same reference!"

    #
    type_dict={}
    for x in bed4_rrna:
        type_dict[x[3]]="rRNA"
    for x in bed4_trna:
        type_dict[x[3]]="tRNA"
    for x in bed4_cds:
        type_dict[x[3]]="CDS"

    bedall=sorted(bed4_rrna+bed4_cds+bed4_trna)

    ####todo: add the position for the non-coding region as NCR, cutoff =gap>20
    ####

    ####

    out_l=[]

    for line in bedall:
        chro, start, end, anno=line
        if type_dict[anno]=="tRNA":

            seq3="tRNA-"+str(SeqUtils.seq3(anno))
            line2w="{start}\t{end}\t{type}\n\t\t\t{key}\t{value}\n".format(
                start=start,end=end, type="tRNA",key="product",value=seq3)

        elif type_dict[anno]=="rRNA":
            line2w="{start}\t{end}\t{type}\n\t\t\t{key}\t{value}\n".format(
                start=start,end=end, type="rRNA",key="product",value=anno)

        elif type_dict[anno]=="CDS":
            line2w_1="{start}\t{end}\t{type}\n\t\t\t{key}\t{value}\n".format(
                start=start,end=end, type="gene",key="gene",value=anno)
            line2w_2="{start}\t{end}\t{type}\n\t\t\t{key1}\t{value1}\n\t\t\t{key2}\t{value2}\n".format(
                start=start,end=end, type="CDS",
                key1="product",value1=anno,
                key2="transl_table",value2=5)
            line2w="".join([line2w_1, line2w_2])

        out_l.append(line2w)

    return out_l


def _cmsearch_wrapper_rrna(fastafile, MITFIPATH=None):

    """
    todo: too slow to be practical, maybe change to INFERNAL 1.1 and try
    mitfi.jar in in $MITFIPATH=./bins
    :return:teh filename of mitfi run
    """
    if MITFIPATH is None:
        path = os.path.dirname(__file__)
        MITFIPATH=os.path.join(path, "bins", "mitfi")
        #print MITFIPATH
    jarfile=os.path.join(MITFIPATH, "mitfi.jar")
    rrna_cm=os.path.join(os.path.dirname(__file__), "bins", "mitfi","r_rna.cm")

    mitfi_cmd = "java -jar {jarfile} -cm {rrna_cm} -top {fastafile}".format(
        jarfile=jarfile, fastafile=fastafile, rrna_cm=rrna_cm)
    rrna_out = myexe(mitfi_cmd)
    print(rrna_out)
    prefix=fastafile.split("/")[-1].split(".")[0]
    with open(prefix+"_rrna.txt", "w") as fw:
        fw.write(rrna_out)

    return prefix+"_rrna.txt"


def get_tnra_pro(fastafile, mitfi_out):
    """
    re_order the mtDNA fastafile, write a new fastefile start with trna P (Phe)

    :param fastafile:
    :param MITFIPATH:
    :return:
    """
    fasta_d=fasta2dic(fastafile)
    if len(fasta_d)!=1:
        print("Check the fasta file and make it continues one!")
        return None

    ###### get the start point of trna_Pro and the strand
    trna_pro_start=None
    trna_pro_socre=0
    strand="+"
    fr=open(mitfi_out,"r")
    for line in fr.readlines():
        terms=line.strip().split("\t")
        header, start, stop, score, evalue, AC, AA, model, strand=terms
        if AA=="P" and float(score)>trna_pro_socre and float(evalue)<=0.001:
            if strand=="+":
                trna_pro_start=int(start)
            elif strand=="-":
                trna_pro_start=int(stop)
            trna_pro_socre=score
    fr.close()

    return (trna_pro_start, strand)


def re_order(fastafile,newstart,strand="+",outfasta=None):
    """
    :param fastafile: the sequence which only have one sequence
    :param newstart: new start point for the fasta sequence, 1 based
    :param strand": "+" or "-"
    :param outfasta:
    :return:
    """
    fasta_d=fasta2dic(fastafile)

    ######
    if len(fasta_d)!=1:
        print("Check the fasta file and make it continues one!")
        return None
    chro=fasta_d.keys()[0]
    seq=fasta_d.values()[0]
    ######

    ###### re-order the new file
    if outfasta is None:
        prefix=fastafile.split("/")[-1].split(".")[0]
        outfasta=prefix+"_ordered.fasta"
    ########

    if strand=="+":
        pass
    if strand=="-":
        newstart=len(seq)-newstart+1
        seq=seq.reverse_complement()
        fasta_d={chro:seq}
        print(newstart)

    frg_1=chr_select(record_dict=fasta_d, chro=chro, start=newstart-1, end=len(seq))[1]
    frg_2=chr_select(record_dict=fasta_d, chro=chro, start=0, end=newstart-1)[1]
    seq_new="".join([frg_1, frg_2])

    with open(outfasta, "w") as fw:
        fw.write(">"+chro+"_re"+"\n")
        fw.write(seq_new)
        fw.write("\n")
    return outfasta


def flow_re_order(fastafile, MITFIPATH=None, outfasta=None):


    mitfi_out=mitfi_wrapper_trna(fastafile, MITFIPATH)
    newstart, strand=get_tnra_pro(fastafile, mitfi_out)
    outfasta=re_order(fastafile, newstart, strand, outfasta)

    return outfasta


#### flow for tbl generation
def pre_fsa(fastafile, spe_name=None, out=None):
    """
    The species name should be provided in the running or
    just rename the fastafile to contain the species name

    prepare the fsa file for the tbl2asn prog
    the fsa file should:
        - have the header line contains the spe information
        - end with ".fsa"
    :param fastafile:
    :param spe_name:
    :param out:
    :return:
    """
    fa_d=fasta2dic(fastafile)
    if spe_name is None:
        spe_name= fastafile.split(".")[0]
    if out is None:
        out = spe_name + ".fsa"
    fastaname=spe_name.replace(" ", "_") # in case the species name contains space

    # give a mtDNA header line for the fsa file, ready for the submission
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


def flow_tbl(fasta_order, cds_ref, r_ref, spe_name=None,MITFIPATH=None):
    #r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
    #cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"

    #fasta_order=flow_re_order(fastafile, MITFIPATH=MITFIPATH)
    prefix=".".join(fasta_order.split(".")[0:-1])

    if spe_name is None:
        spe_name=prefix

    trnafile=mitfi_wrapper_trna(fasta_order, MITFIPATH, prefix=prefix)
    bed4_trna=trnafile_parser(trnafile)
    bed4_rrna=rrna_tbl_get(fasta_order, r_ref)
    bed4_cds=cds_tbl_get(fasta_order, cds_ref)
    tbl_out=tbl_format(bed4_rrna, bed4_cds, bed4_trna)


    fsa_name=pre_fsa(fasta_order, spe_name)

    headstr=">Feature\t"+fsa_name+"\n"
    tbl_out=[headstr]+tbl_out
    tbl_name=spe_name+".tbl"

    with open(tbl_name, "w") as fw:
        fw.write("".join(tbl_out))

    print("Using {fsa_name} and {tbl_name} for your tbl2asn submission.".format(
        fsa_name=fsa_name, tbl_name=tbl_name
    ))

    return fsa_name, tbl_name

if __name__=="__main__":
    os.chdir("/home/zhaolab1/data/mitosra/dna/wkdir/1561998")
    r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
    cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"
    out=flow_re_order("1561998.fasta")
    flow_tbl(out, cds_ref, r_ref, spe_name="Caenorhabditis tropicalis")
