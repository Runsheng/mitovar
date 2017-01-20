#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2017/1/19 10:33
# @Author  : Runsheng     
# @File    : anno_test.py.py

from anno import *
from glob import glob

def flow_tbl():
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

if __name__=="__main__":
    print flow_tbl()
    pre_submission()
    print "done"
