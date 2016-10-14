#!/usr/bin/env python
#
# File: rarepedsim_utilities.py
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zong-Xiao He"
__date__ = "11/30/15"

# from pedigree import *
from catt_weight import *
from disequilibrium_test import *
import operator


def read_pedfile(pedfile):
    fam2ped = {}
    sid2geno = {}
    with open(pedfile) as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            items = line.strip().split()
            fid = items[0]
            info = items[1:6]
            if fid not in fam2ped:
                fam2ped[fid] = [info]
            else:
                fam2ped[fid].append(info)

            sid = items[1]
            # output of raresimped: 1 for wildtype, 2 for minor, and two columns for a genotype
            raw_geno = [int(i) for i in items[6:]]
            geno = []
            allmissing = True
            for i,j in zip(raw_geno[0::2], raw_geno[1::2]):
                if i == 0 or j == 0:
                    geno.append(-9)
                else:
                    allmissing = False
                    geno.append(i + j - 2)

            if not allmissing:
                sid2geno[sid] = geno
            # geno = [int(i) fo r i in items[6:]]
            # geno = [int(i) - 1 for i in items[6:]]
            # geno = [i+j for i,j in zip(geno[0::2], geno[1::2])]

    peds = []
    for fam in fam2ped:
        ped = Pedigree(fam, fam2ped[fam])
        peds.append(ped)
    return peds, sid2geno


def count_population_ctrl(pedfile):
    minor, count = None, None
    with open(pedfile) as f:
        for i, line in enumerate(f.readlines()):
            items = line.strip().split()
            geno = [int(i) - 1 for i in items[6:]]
            geno = [i+j for i,j in zip(geno[0::2], geno[1::2])]
            m = [i if i >= 0 else 0 for i in geno]
            c = [2 if i >= 0 else 0 for i in geno]
            minor = m if minor is None else map(operator.add, minor, m)
            count = c if count is None else map(operator.add, count, c)
    return [[i, j] for (i,j) in zip(minor, count)]


def count_founder(peds, sid2geno):
    num_sites = len(sid2geno.values()[0])
    minor, count = [0]*num_sites, [0]*num_sites
    for ped in peds:
        for founder in ped.founders:
            if founder not in sid2geno:
                continue

            for i, g in enumerate(sid2geno[founder]):
                if g < 0:
                    continue
                minor[i] += g
                count[i] += 2
    return [[m, c] for m,c in zip(minor, count)]


def popctrl_weight(founder_count, ctrl_count):
    weights = []
    for case, ctrl in zip(founder_count, ctrl_count):
        weights.append(catt_weight(case, ctrl))
    return weights


def read_mafs(sfsfile):
    mafs = []
    with open(sfsfile) as f:
        for line in f.readlines():
            if line.startswith("#"):
                continue
            items = line.strip().split()
            mafs.append(float(items[3]))
    return mafs


def rarepedsim_2_rvdt(case_pedfile, sfsfile):
    peds, sid2geno = read_pedfile(case_pedfile)

    # if ctrl_pedfile is None:
    #     weights = [1.0] * len(sid2geno.values()[0])
    # else:
    #     founder_count = count_founder(peds, sid2geno)
    #     popctrl_count = count_population_ctrl(ctrl_pedfile)
    #     weights = popctrl_weight(founder_count, popctrl_count)

    mafs = read_mafs(sfsfile)
    # if len(mafs) != len(weights):
    #     raise ValueError("Mafs length, weights length, and genotype length don't match")
    return (peds, sid2geno, mafs)