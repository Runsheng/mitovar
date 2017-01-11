#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/8 16:58
# @Author  : Runsheng     
# @File    : utils.py

import subprocess
import sys
import signal
import os
import fnmatch

from Bio import SeqIO

def myexe(cmd, timeout=0):
    """
    a simple wrap of the shell
    mainly used to run the bwa mem mapping and samtool orders
    """
    def setupAlarm():
        signal.signal(signal.SIGALRM, alarmHandler)
        signal.alarm(timeout)

    def alarmHandler(signum, frame):
        sys.exit(1)

    proc=subprocess.Popen(cmd, shell=True, preexec_fn=setupAlarm,
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=os.getcwd())
    out, err=proc.communicate()
    print err
    return out, err, proc.returncode


def fasta2dic(fastafile):
    """
    Give a fasta file name, return a dict contains the name and seq
    Require Biopython SeqIO medule to parse the sequence into dict, a large genome may take a lot of RAM
    """
    if ".gz" in fastafile:
        handle=gzip.open(fastafile, "rU")
    else:
        handle=open(fastafile, "rU")
    record_dict=SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
    handle.close()
    return record_dict


def chr_select(record_dict, chro, start,end):
    """
    Note the start and end is 0 based
    give the name of refdic, and the chr, start and end to be used
    return the name and sequence (both as str)
    for example, chrcut(record_dict, "I", 100,109) returns
     ("I:100_109","AAAAAAAAAA")
    """
    name=chro+ ":"+str(start)+"_"+str(end)
    seq=str(record_dict[chro][start:end].seq)
    return name,seq


def dic2fasta(record_dict,out="record_dict.fasta"):
    """
    Write a record_dict of SeqIO sequence dict back to fasta file
    :param record_dict:
    :param out:
    :return:
    """
    with open(out,"w") as f:
        for record in sorted(record_dict.keys()):
            name=record
            seq=str(record_dict[name].seq)
            f.write(">")
            f.write(name)
            f.write("\n")
            f.write(seq)
            f.write("\n")


def myglob(seqdir, word):
    """
     to write a glob for python2 for res-glob
    """
    matches=[]
    for root, dirnames, filenames in os.walk(seqdir):
         for filename in fnmatch.filter(filenames, word):
            matches.append(os.path.join(root, filename))
    return matches