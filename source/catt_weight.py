#!/usr/bin/env python
#
# File: trend_test.py
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zong-Xiao He"
__date__ = "11/17/15"

import numpy as np

def cochran_armitage_trend_test(row1, row2):
    """
    Cochran-Armitage Trend Test statistics for 2xc tables

    Reference: https://code.google.com/p/glu-genetics/source/browse/glu/lib/ca.py?r=9947b985c02808b240f6a706bbe6f4af025a757f
    """
    x = np.asarray([row1, row2], dtype=int)
    if len(x.shape) != 2 or x.shape[0] != 2:
        raise ValueError,'cochran_armitage_trend_test requires a 2xc table'

    n_i   = x.sum(axis=0)
    n     = n_i.sum()
    r_i   = np.arange(len(n_i))
    r_bar = float((n_i*r_i).sum())/n
    s2    = (n_i*(r_i-r_bar)**2).sum()
    p1    = float(x[0].sum())/n

    t = (x[0]*(r_i-r_bar)/(p1*(1-p1)*s2)**0.5).sum()

    return t


def catt_weight(founder, pop_ctrl):
    if founder[1] == 0 or pop_ctrl[1] == 0:
        raise ValueError, 'can not calculate weight for uninformative variant site'

    if founder[0] == 0 and pop_ctrl[0] == 0:
        return 0
    else:
        catt = cochran_armitage_trend_test(founder, pop_ctrl)
        dir = 1 if founder[0]*1.0/founder[1] > pop_ctrl[0]*1.0/pop_ctrl[1] else -1
        return abs(catt) * dir


if __name__ == '__main__':
    # should return -4.7918
    print cochran_armitage_trend_test([26,26,23,18,9],[6,7,9,14,23])

    print catt_weight([2, 12], [200, 2000])
    print catt_weight([0,32],[1,32])
    print catt_weight([0, 6000], [1, 10000])