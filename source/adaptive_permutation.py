#!/usr/bin/env python
#
# File: adaptive_permutation.py
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zong-Xiao He"
__date__ = "11/5/15"

import random
import math
from scipy.stats import norm


class AdaptivePermutation(object):
    def __init__(self, test_score, adapt_iter=1000, max_iter=2e6, alpha=2.5e-6):
        self.test_score = test_score
        self.num_test = len(self.test_score)
        self.adapt_iter = adapt_iter
        self.max_iter = max_iter
        self.alpha = alpha
        self.count = [0] * self.num_test
        self.pvalue = [0] * self.num_test
        self.iter = 0

    def get_pvalue(self):
        return self.pvalue

    def update_count(self, permutation_score):
        self.iter += 1
        for i in range(self.num_test):
            diff = permutation_score[i] - self.test_score[i]
            if diff > 1e-6:
                self.count[i] += 1
            elif -1e-6 <= diff <= 1e-6:
                self.count[i] += random.randint(0, 1)
            else:
                continue

        # print permutation_score, self.test_score, "-->", self.count
        should_break = False
        if self.iter % self.adapt_iter == 0 or self.iter >= self.max_iter:
            should_break = self.__update_pvalue()

        if should_break is True or self.iter >= self.max_iter:
            return True
        else:
            return False

    def __update_pvalue(self):
        should_break = True
        for i, c in enumerate(self.count):
            if c == self.iter:
                self.pvalue[i] = 1
                continue
            # test the hypothesis: pval is close to alpha
            pval = (c + 1.0) / (self.iter + 1.0)
            sigma = math.sqrt(pval * (1 - pval) / self.iter)
            cvalue = norm.ppf(0.975, scale=sigma)

            self.pvalue[i] = pval
            if (pval - self.alpha) < cvalue:
                should_break = False

        return should_break


if __name__ == '__main__':
    permut = AdaptivePermutation([3,4], 200, 200000, 0.00003)
    count = 0
    while True:
        count += 1
        score = random.gauss(0, 1)
        sbreak = permut.update_count([score, score])
        if sbreak:
            print count, permut.get_pvalue()
            break