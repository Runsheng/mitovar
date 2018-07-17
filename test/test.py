#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 6/27/2018 2:12 PM
# @Author  : Runsheng     
# @File    : test.py

"""
The test for the mitovar main function, currently:
1 assemble
2 anno
"""
from utils import myexe
import os


def test_anno():
    """
    mimic the cmd usage of mitovar
    :return:
    """
    path = os.path.dirname(__file__)
    print("Anno test may take around 5 min, please be patient.")
    cmd="""
    cd {path}
    ../mitovar.py anno -f 1561998.fasta -c celmt_p.fasta -r rrna.fasta -s ctro
    """
    out=myexe(cmd)

    if "ctro.tbl" in out:
        print("anno test pass")
    else:
        print("Error, please chech the log")

def test_assemble():
    pass



if __name__=="__main__":
    test_anno()