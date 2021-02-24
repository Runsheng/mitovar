#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/1/19 10:33
# @Author  : Runsheng     
# @File    : anno_test.py.py

from anno import *
from glob import glob

def flow_tbl(work_dir, r_ref, cds_ref, MITHIPATH=None):
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/trna"
    os.chdir(wkdir)
    r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
    cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"

    trna_files=sorted(glob("*_trna.txt"))
    fastafiles=sorted(glob("*_ordered.fasta"))

    for seqfile, trnafile in zip(fastafiles,trna_files):
        prefix=seqfile.split(".")[0]
        bed4_rrna=rrna_tbl_get(seqfile, r_ref)
        bed4_cds=cds_tbl_get(seqfile, cds_ref)
        bed4_trna=trnafile_parser(trnafile)
        tbl_out=tbl_format(bed4_rrna, bed4_cds, bed4_trna)
        pre_fsa(seqfile)
        try:
            taxid = int(prefix.split("_")[0])
        except ValueError:
            spe_name=prefix

        spe_name = taxid_d[taxid]
        fastaname = spe_name.replace(" ", "_")
        headstr=">Feature\t"+fastaname+"\n"
        tbl_out=[headstr]+tbl_out

        with open(prefix+".tbl", "w") as fw:
            fw.write("".join(tbl_out))

def flow_tbl_one(fastafile, r_ref, cds_ref, work_dir=None, MITHIPATH=None):

    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/trna"
    os.chdir(wkdir)
    r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
    cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"

    trna_files=sorted(glob("*_trna.txt"))
    fastafiles=sorted(glob("*_ordered.fasta"))

    for seqfile, trnafile in zip(fastafiles,trna_files):
        prefix=seqfile.split(".")[0]
        bed4_rrna=rrna_tbl_get(seqfile, r_ref)
        bed4_cds=cds_tbl_get(seqfile, cds_ref)
        bed4_trna=trnafile_parser(trnafile)
        tbl_out=tbl_format(bed4_rrna, bed4_cds, bed4_trna)
        pre_fsa(seqfile)
        try:
            taxid = int(prefix.split("_")[0])
        except ValueError:
            spe_name=prefix

        spe_name = taxid_d[taxid]
        fastaname = spe_name.replace(" ", "_")
        headstr=">Feature\t"+fastaname+"\n"
        tbl_out=[headstr]+tbl_out

        with open(prefix+".tbl", "w") as fw:
            fw.write("".join(tbl_out))


def pre_submission():
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/trna/"
    os.chdir(wkdir)
    fastafiles =sorted(glob("*_ordered.fsa"))
    for i in fastafiles:
        cmd="cp {i} ./submission/{i}".format(i=i)
        myexe(cmd)
    myexe("cp *.tbl ./submission/")
    myexe("tbl2asn -p ./submission -t ./submission/author4.sbt" )

taxid_d={31234: 'Caenorhabditis remanei',
 135651: 'Caenorhabditis brenneri',
 281681: 'Caenorhabditis plicata',
 281687: 'Caenorhabditis japonica',
 497829: 'Caenorhabditis sinica',
 860376: 'Caenorhabditis angaria',
 1094320: 'Caenorhabditis sp. 1 KK-2011',
 1094321: 'Caenorhabditis doughertyi',
 1094323: 'Caenorhabditis virilis',
 1094326: 'Caenorhabditis wallacei',
 1094327: 'Caenorhabditis nouraguensis',
 1094328: 'Caenorhabditis macrosperma',
 1094331: 'Caenorhabditis guadeloupensis',
 1094335: 'Caenorhabditis afra',
 1561998: 'Caenorhabditis tropicalis',
 1630362: 'Caenorhabditis castelli',
 1729975: 'Caenorhabditis sp. 38 MB-2015'}


def pre_fsa(fastafile):
    fa_d=fasta2dic(fastafile)
    prefix = fastafile.split(".")[0]
    out = prefix + ".fsa"
    taxid=int(prefix.split("_")[0])
    spe_name=taxid_d[taxid]
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


def test_flow_one(fasta):
    """
    test a single one, from
    :return:
    """
    pass


def test_flowexon():
    from glob import glob
    import os
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/exon"
    os.chdir(wkdir)
    target = "./ref/celmt_p.fasta"
    query_list = glob("*_s.fasta")
    print(query_list)
    for query in query_list:
        print query
        flow_exon(query, target)
    myexe("""
     cat *_p_corr.fasta> m20_p.fasta
     cat *_cds_corr.fasta >m20_cds.fasta
     """)


def test_corr():
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/exon/round0"
    os.chdir(wkdir)
    _exon_corrector("281687_s_cds.fa", "281687_s_p.fa")

    # test_flowexon()


def test_cog():
    get_cogfile("m20_p.fasta")


def test_genbank_parser():
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/exon"
    os.chdir(wkdir)


def test_miftf():
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/exon"
    os.chdir(wkdir)
    print mitfi_wrapper_trna("/home/zhaolab1/myapp/mitovar/bins/mitfi/test/celmt_re.fasta")


def test_reorder():
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/exon"
    os.chdir(wkdir)
    re_order("/home/zhaolab1/myapp/mitovar/bins/mitfi/test/celmt_re.fasta", 13794, "-")
    # the re_order need to be same as org fasta file


def test_getpro():
    print get_tnra_pro("/home/zhaolab1/myapp/mitovar/bins/mitfi/test/celmt_re.fasta",
                       "/home/zhaolab1/myapp/mitovar/celmt_re_trna.txt")


def test_re_flow():
    from pprint import pprint
    wkdir = "/home/zhaolab1/data/mitosra/dna/anno/trna"
    os.chdir(wkdir)
    from glob import glob
    from utils import parmap
    # files=glob("*_s.fasta")
    # parmap(re_order_flow, files, 32)
    files = glob("*_ordered.fasta")
    # parmap(mitfi_wrapper_trna, files, 32)
    r_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/rrna.fasta"
    cds_ref = "/home/zhaolab1/data/mitosra/dna/anno/exon/ref/celmt_p.fasta"
    for i in files:
        pprint(rrna_tbl_get(i, r_ref))
        pprint(cds_tbl_get(i, cds_ref))

if __name__=="__main__":
    print flow_tbl()
    pre_submission()
    print "done"
