#!/usr/bin/env python
#
# File: rvdt_utilities.py
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zong-Xiao He"
__date__ = "11/30/15"

# from pedigree import *
from disequilibrium_test import *
import sys

def read_weights(wfile):
    weights = []
    with open(wfile) as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue

            try:
                weights.append(float(line.strip()))
            except:
                sys.exit("Weights needs to be float number")
    return weights


def read_genofile(genofile):
    sid2geno = {}
    with open(genofile) as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            items = line.strip().split()
            sid = items[0]
            geno = [int(i) for i in items[1:]]
            sid2geno[sid] = geno

    return sid2geno


def read_pedfile(pedfile):
    fam2ped = {}
    with open(pedfile) as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            items = line.strip().split()
            # print items
            fid = items[0]
            info = items[1:]
            if fid not in fam2ped:
                fam2ped[fid] = [info]
            else:
                fam2ped[fid].append(info)

    peds = []
    for fam in fam2ped:
        ped = Pedigree(fam, fam2ped[fam])
        peds.append(ped)
    return peds

