#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2016/12/8 16:58
# @Author  : Runsheng     
# @File    : utils.py

import subprocess
import sys
import signal
import os


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
