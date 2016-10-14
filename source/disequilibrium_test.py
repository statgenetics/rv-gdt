#!/usr/bin/env python
#
# File: disequilibrium_test.py
#
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zong-Xiao He"
__date__ = "11/2/15"

# import random
import operator
import numpy as np
# import math
from adaptive_permutation import *
# from catt_weight import *
import scipy.stats as stat
import json
from pedigree import *
# from rarepedsim_utilities import *
import statsmodels.api as sm

# Genotype coding
# -9 missing
#  0 wildtype
#  1 hetero
#  2 homo


class DisequilibriumTest(object):
    """Rare Variant Extesnion of Disequilibrium Test 

    Attributes:
        genotype data: subjet id -> genotype
        family structure: a list of Pedigree object
    """

    def __init__(self, peds, sid2geno, mafs=None, weights=None):
        self.pedigree = peds  # list of Pedigree object
        # self.num_sites = len(mafs)
        # self.weights = weights if weights is not None else [1]*self.num_sites
        # self.pop_mafs = mafs

        # self.founder_count = []
        self._pre_process(sid2geno, mafs, weights)
        # self._popctrl_weight()

        # test options
        self.test_options = {}
        self.tests = []
        self.vt_mafs = None

    def remove_founders(self):
        for ped in self.pedigree:
            for founder in ped.founders:
                if founder in self.sid2geno:
                    self.sid2geno.pop(founder, None)
                    # del self.sid2geno[founder]
        return

    def _pre_process(self, sid2geno, mafs, weights):
        """
        remove uninformative variant sites
        Updated: remove this function, in some case, there is no founder genotypes
        """
        # TODO: needs to handle data without founder
        self.uninformative_data = False
        num_sites = len(sid2geno.values()[0])
        minor, total = [0]*num_sites, [0]*num_sites
        for ped in self.pedigree:
            for founder in ped.founders:
                if founder not in sid2geno:
                    continue

                for i, g in enumerate(sid2geno[founder]):
                    if g < 0:
                        continue
                    minor[i] += g
                    total[i] += 2

        # informative_site = []
        mafs_in_founder = []
        for m,t in zip(minor, total):
            if m == 0 or t == 0:
                # informative_site.append(False)
                mafs_in_founder.append(0)
            else:
                # informative_site.append(True)
                mafs_in_founder.append(m * 1.0 / t)

        # use maf in founder is pop maf is not given
        self.pop_mafs = mafs_in_founder if (mafs is None) else mafs
        # use uniform weights if weight is not given
        self.weights = [1]*num_sites if (weights is None) else weights

        self.sid2geno = sid2geno

        return

        # # remove variant sites with weight = 0
        # for i, w in enumerate(my_weights):
        #     if w == 0:
        #         informative_site[i] = False
        #
        # self.pop_mafs = [m for m,i in zip(my_mafs, informative_site) if i]
        # self.weights = [w for w,i in zip(my_weights, informative_site) if i]
        #
        # if len(self.pop_mafs) == 0:
        #     self.uninformative_data = True
        #     return
        # else:
        #     self.uninformative_data = False
        #
        # self.sid2geno = {}
        # for sid in sid2geno:
        #     geno = [g for informative, g in zip(informative_site, sid2geno[sid]) if informative]
        #     self.sid2geno[sid] = geno

    def set_test_options(self, args):

        # args.maf_cutoff, args.unweighted, args.use_vt, args.analytical, args.vpdt
        # args.adapt, args.max_iter, args.alpha
        self.test_options["maf_cutoff"] = args.maf_cutoff
        self.test_options["include_unweighted"] = args.unweighted
        self.test_options["include_vt"] = args.use_vt
        self.test_options["include_analytical_pvalues"] = args.analytical
        self.test_options["adapt_iter"] = args.adapt_iter
        self.test_options["max_iter"] = args.max_iter
        self.test_options["alpha"] = args.alpha
        self.test_options["vpdt"] = args.vpdt
        # self.test_options["shuffle_method"] = shuffle_method

        # all tests
        gdt = ["RVGDT"]
        if args.use_vt:
            gdt.append("RVGDT-VT")
        if args.unweighted:
            gdt += [i + "-UnWeighted" for i in gdt]

        if args.vpdt:
            pdt = ["RVPDT"]
            if args.use_vt:
                pdt.append("RVPDT-VT")
            if args.unweighted:
                pdt += [i + "-UnWeighted" for i in pdt]
            self.tests = pdt + gdt
        else:
            self.tests = gdt
        return


    def analyze(self):
        if self.test_options["vpdt"]:
            tests, pvals = self.rvpdt()
        else:
            tests, pvals = self.rvgdt()

        if pvals is  None:
            pvals = ['none' for t in tests]

        return tests, pvals


    def rvpdt(self):
        """
        Analyze the family data with RVPDT, return name of tests and corresponding pvalues
        """
        # if self.uninformative_data:
        #     return None, None

        # Updated: Aug 15 2016
        # rvpdt can not handle missing family member correctly, if there is missing 'founder' within the nuclear family
        # we need to exclude this nuclear family
        genotyped = self.sid2geno.keys()
        for ped in self.pedigree:
            ped.preprocess_rvpdt(genotyped)

        tests = ["RVPDT"]
        if self.test_options["include_vt"]:
            tests.append("RVPDT-VT")
        if self.test_options["include_unweighted"]:
            tests += [i + "-UnWeighted" for i in tests]

        # original statistics
        pdt_Sij = self._rvpdt_statistic_matrix()

        # we need to calculate the mafs used in vt approach
        if self.test_options["include_vt"] is True:
            self.vt_mafs = self._maf_in_affected(self.sid2geno)
            scores = self.__weight_and_collapse_vt(pdt_Sij)
        else:
            scores = self.__weight_and_collapse(pdt_Sij)
        # print "\t".join([str(s) for s in scores])

        if scores is None:
            return tests, None

        # analytical solution if max_iter is 0
        if self.test_options["max_iter"] == 0:
            alltests = ["Analytical-" + i for i in self.tests]
            pvalues = [stat.norm.sf(s) for s in scores]
            return alltests, pvalues

        permutation = AdaptivePermutation(scores, self.test_options["adapt_iter"],
                                          self.test_options["max_iter"],
                                          self.test_options["alpha"])
        while True:
            self._genotype_shuffle()
            # vt_mafs = self._maf_in_affected(sid2geno)
            pdt_Sij = self._rvpdt_statistic_matrix()

            if self.test_options["include_vt"] is True:
                pscores = self.__weight_and_collapse_vt(pdt_Sij)
            else:
                pscores = self.__weight_and_collapse(pdt_Sij)

            should_break = permutation.update_count(pscores)
            if should_break:
                break

        if self.test_options["include_analytical_pvalues"]:
            alltests = ["Analytical-" + i for i in tests] + ["Empirical-" + i for i in tests]
            pvalues = [stat.norm.sf(s) for s in scores] + permutation.get_pvalue()
            return alltests, pvalues
        else:
            return ["Empirical-" + i for i in tests], permutation.get_pvalue()


    def rvgdt(self):
        """
        Analyze the family data with RVGDT, return name of tests and corresponding pvalues
        """
        # if self.uninformative_data:
        #     return None, None

        tests = ["RVGDT"]
        if self.test_options["include_vt"]:
            tests.append("RVGDT-VT")
        if self.test_options["include_unweighted"]:
            tests += [i + "-UnWeighted" for i in tests]

        # original statistics
        gdt_Sij = self._rvgdt_statistic_matrix()

        # we need to calculate the mafs used in vt approach
        if self.test_options["include_vt"] is True:
            self.vt_mafs = self._maf_in_affected(self.sid2geno)
            scores = self.__weight_and_collapse_vt(gdt_Sij)
        else:
            scores = self.__weight_and_collapse( gdt_Sij)
        # print "\t".join([str(s) for s in scores])

        if scores is None:
            return tests, None

        # analytical solution if max_iter is 0
        if self.test_options["max_iter"] == 0:
            alltests = ["Analytical-" + i for i in self.tests]
            pvalues = [stat.norm.sf(s) for s in scores]
            return alltests, pvalues

        permutation = AdaptivePermutation(scores, self.test_options["adapt_iter"],
                                          self.test_options["max_iter"],
                                          self.test_options["alpha"])
        while True:
            self._phenotype_shuffle()
            # vt_mafs = self._maf_in_affected(sid2geno)
            gdt_Sij = self._rvgdt_statistic_matrix()

            if self.test_options["include_vt"] is True:
                pscores = self.__weight_and_collapse_vt(gdt_Sij)
            else:
                pscores = self.__weight_and_collapse(gdt_Sij)

            should_break = permutation.update_count(pscores)
            if should_break:
                break

        if self.test_options["include_analytical_pvalues"]:
            alltests = ["Analytical-" + i for i in tests] + ["Empirical-" + i for i in tests]
            pvalues = [stat.norm.sf(s) for s in scores] + permutation.get_pvalue()
            return alltests, pvalues
        else:
            return ["Empirical-" + i for i in tests], permutation.get_pvalue()

    # def compare_nosibpair(self, maf_cutoffs=[0.01], adapt_iter=1000, max_iter=2e6, alpha=2.5e-6,
    #                       include_noweight=True, include_analytical=True):
    #     if self.uninformative_data:
    #         return None, None
    #
    #     # original statistics
    #     scores = []
    #     # vt_mafs = self._maf_in_affected(self.sid2geno)
    #     Sij_sib, Sij_nosib = self._calculate_statistic_matrix2(self.sid2geno)
    #
    #     for cutoff in maf_cutoffs:
    #         scores += self.__weight_and_collapse(Sij_sib,  cutoff, include_noweight)
    #         scores += self.__weight_and_collapse(Sij_nosib,  cutoff, include_noweight)
    #     # print "\t".join([str(s) for s in scores])
    #
    #     permutation = AdaptivePermutation(scores, adapt_iter, max_iter, alpha)
    #     while True:
    #         # shuffle
    #         sid2geno = self._shuffle()
    #         # vt_mafs = self._maf_in_affected(sid2geno)
    #         Sij_sib, Sij_nosib = self._calculate_statistic_matrix2(sid2geno)
    #
    #         pscores = []
    #         for cutoff in maf_cutoffs:
    #             pscores += self.__weight_and_collapse(Sij_sib,  cutoff, include_noweight)
    #             pscores += self.__weight_and_collapse(Sij_nosib,  cutoff, include_noweight)
    #
    #         # print "\t".join([str(s) for s in pscores])
    #         should_break = permutation.update_count(pscores)
    #         if should_break:
    #             break
    #
    #     # organize output: return analytical p values and empirical p values
    #     tests = []
    #     for c in maf_cutoffs:
    #         tests.append("{0}-T{1}".format(self.test_name, str(c)[2:]))
    #         if include_noweight:
    #             tests.append("{0}-NoWeight-T{1}".format(self.test_name, str(c)[2:]))
    #
    #         tests.append("{0}-T{1}-NoSib".format(self.test_name, str(c)[2:]))
    #         if include_noweight:
    #             tests.append("{0}-NoWeight-T{1}-NoSib".format(self.test_name, str(c)[2:]))
    #
    #     if include_analytical:
    #         alltests = ["Analytical-" + i for i in tests] + ["Empirical-" + i for i in tests]
    #         pvalues = [stat.norm.sf(s) for s in scores] + permutation.get_pvalue()
    #         return alltests, pvalues
    #     else:
    #         return ["Empirical-" + i for i in tests], permutation.get_pvalue()

    def __weight_and_collapse(self, pdt_matrix):
        # uninformative
        if len(pdt_matrix) == 0:
            return None

        # print "Matrix:", matrix
        should_analyze = [True if (m <= self.test_options["maf_cutoff"]) else False for m in self.pop_mafs]

        # pdt
        weighted_matrix = [[s*w for s, w in zip(row, self.weights)] for row in pdt_matrix]
        Di = [sum([s if a else 0 for s, a in zip(row, should_analyze)]) for row in weighted_matrix]
        weighted_score = self.__test_statistic(Di)

        if self.test_options["include_unweighted"]:
            Di = [sum([s if a else 0 for s, a in zip(row, should_analyze)]) for row in pdt_matrix]
            unweighted_score = self.__test_statistic(Di)
            dt = [weighted_score, unweighted_score]
        else:
            dt = [weighted_score]

        return dt

    # similar to pdt and gdt
    def __weight_and_collapse_vt(self, pdt_matrix):
        # uninformative
        if len(pdt_matrix) == 0:
            return None

        # print "Matrix:", matrix
        sorted_mafs = []
        for (m,n) in zip(self.pop_mafs, self.vt_mafs):
            if m <= self.test_options["maf_cutoff"] and n > 0:
                sorted_mafs.append(n)
        sorted_mafs = list(set(sorted_mafs))
        sorted_mafs.sort()

        # weighted
        weighted_pdt = [[s*w for s, w in zip(row, self.weights)] for row in pdt_matrix]

        pdt_vt, pdt_novt = float('-inf'), float('-inf')
        pdt_vt_unweighted, pdt_novt_unweighted = float('-inf'), float('-inf')

        for cutoff in sorted_mafs:
            should_analyze = [True if (m <= self.test_options["maf_cutoff"] and n <= cutoff) else False for (m,n) in zip(self.pop_mafs, self.vt_mafs)]

            # pdt
            Di = [sum([s if a else 0 for s, a in zip(row, should_analyze)]) for row in weighted_pdt]
            score = self.__test_statistic(Di)
            pdt_vt = max(score, pdt_vt) # max of variable threshold
            pdt_novt = score # including all should analyzed variant = no vt score

            if self.test_options["include_unweighted"]:
                Di = [sum([s if a else 0 for s, a in zip(row, should_analyze)]) for row in pdt_matrix]
                unweighted_score = self.__test_statistic(Di)
                pdt_vt_unweighted = max(unweighted_score, pdt_vt_unweighted) # max of variable threshold
                pdt_novt_unweighted = unweighted_score # including all should analyzed variant = no vt score


        if self.test_options["include_unweighted"]:
            return [pdt_novt, pdt_vt, pdt_novt_unweighted, pdt_vt_unweighted]
        else:
            return [pdt_novt, pdt_vt]


    @staticmethod
    def __test_statistic(score):
        if sum(score) == 0:
            # TODO:
            return 0
        else:
            return sum(score) * 1.0 / math.sqrt(sum([i*i for i in score]))

    def _maf_in_affected(self, sid2geno):
        """
        MAF in affected subjects, will be used in Variable Threshold
        """
        not_missing_count, allele_count = None, None
        for ped in self.pedigree:
            for sid in ped.affected:
                if sid in sid2geno:
                    geno  = [i if i >= 0 else 0 for i in sid2geno[sid]]
                    count = [1 if i >= 0 else 0 for i in sid2geno[sid]]
                    allele_count = geno if allele_count is None else map(operator.add,
                                                                         allele_count,
                                                                         geno)
                    not_missing_count = count if not_missing_count is None else map(operator.add,
                                                                                    not_missing_count,
                                                                                    count)
        return [i * 0.5 / j if j > 0 else 0 for (i,j) in zip(allele_count, not_missing_count)]


    def _phenotype_shuffle(self):
        """
        RV-GDT only: Shuffle the phenotype
        """
        for ped in self.pedigree:
            ped.shuffle_phenotype()
        return

    def _genotype_shuffle(self):
        """
        RV-PDT only: Shuffle the genotype of non-founder, by randomly mating the parents
        """
        shuffled_sid2geno = {}
        for ped in self.pedigree:
            shuffled = []
            for s in ped.founders:
                shuffled.append(s)
                if s in self.sid2geno:
                    shuffled_sid2geno[s] = self.sid2geno[s]
                else:
                    # If one founder doesn't have genotype, we should generate the genotype from population maf.
                    # updated: Aug 19, 2016 [the following command will not be executed
                    shuffled_sid2geno[s] = self.__random_genotype()
                    pass

            for s in ped.non_founders:
                if s not in shuffled:
                    self._shuffle_subject(ped.sid2subject[s], shuffled, shuffled_sid2geno)

        # we need to remove missing subjects in shuffled_sid2geno
        for sid in shuffled_sid2geno:
            if sid in self.sid2geno:
                self.sid2geno[sid] = shuffled_sid2geno[sid]
        return
        # for s in shuffled_sid2geno:
        #     print s, self.sid2geno[s], shuffled_sid2geno[s]
        # return shuffled_sid2geno

    def __random_genotype(self):
        geno = []
        for m in self.pop_mafs:
            g = 0
            if random.random() < m:
                g += 1
            if random.random() < m:
                g += 1
            geno.append(g)
        return geno

    def _shuffle_subject(self, sub, shuffled, shuffled_sid2geno):
        fat, mot = sub.father, sub.mother
        if fat.sid not in shuffled:
            self._shuffle_subject(fat, shuffled, shuffled_sid2geno)

        if mot.sid not in shuffled:
            self._shuffle_subject(mot, shuffled, shuffled_sid2geno)

        fgeno = shuffled_sid2geno[fat.sid]
        mgeno = shuffled_sid2geno[mot.sid]
        geno = []
        for i,j in zip(fgeno, mgeno):
            if i == -9 or j == -9:
                geno.append(-9)
                continue

            g = 0
            if i == 2:
                g += 1
            elif i == 1:
                g += random.randint(0, 1)

            if j == 2:
                g += 1
            elif j == 1:
                g += random.randint(0, 1)
            geno.append(g)

        shuffled_sid2geno[sub.sid] = geno
        shuffled.append(sub.sid)
        return

    def _rvpdt_statistic_matrix(self):
        """
        Calculate PDT socre matrix: Dij (i for pedigree, j for variant)
        Updated: 2015-10-25 remove (n_T + n_S) in denominator
        """
        Dij = []
        for ped in self.pedigree:
            informative_pedigree = False
            Xt, Xs = [0]*len(self.pop_mafs), [0]*len(self.pop_mafs)
            # Nt, Ns = [0]*len(self.pop_mafs), [0]*len(self.pop_mafs)
            # #(minor allele transmitted) - #(minor allele not transmitted) from the heterozygous
            for fat, mot, kid in ped.get_triads():
                if self.sid2geno.get(fat) is None or self.sid2geno.get(mot) is None or self.sid2geno.get(kid) is None:
                    continue

                informative_pedigree = True
                # print "---\nTrio", fat, mot, kid
                xt, nt = self.__transmission_count(self.sid2geno[fat], self.sid2geno[mot], self.sid2geno[kid])
                # print Xt, xt
                Xt = map(operator.add, Xt, xt)
                # Nt = map(operator.add, Nt, nt)
                # print Xt

            # #(minor allele in affected sib) - #(minor allele in unaffected sib)
            for aff, unaff in ped.get_DSPs():
                if self.sid2geno.get(aff) is None or self.sid2geno.get(unaff) is None:
                    continue

                informative_pedigree = True
                # print "---\nSib", aff, unaff
                xs, ns = self.__discordant_count(self.sid2geno[aff], self.sid2geno[unaff])
                Xs = map(operator.add, Xs, xs)
                # Ns = map(operator.add, Ns, ns)
                # print Xs
            # print Xs, Ns

            if informative_pedigree:
                dij = [(xt + xs) * 1.0 for (xt, xs) in zip(Xt, Xs)]
                Dij.append(dij)
        return Dij

    # def _rvpdt_statistic_matrix2(self):
    #     """
    #     Calculate Dij with / without sibpair
    #     """
    #     Dij, nosib_Dij = [], []
    #     for ped in self.pedigree:
    #         Xt, Xs = [0]*len(self.pop_mafs), [0]*len(self.pop_mafs)
    #         # Nt, Ns = [0]*len(self.pop_mafs), [0]*len(self.pop_mafs)
    #         # #(minor allele transmitted) - #(minor allele not transmitted) from the heterozygous
    #         for fat, mot, kid in ped.get_all_triad():
    #             if self.sid2geno.get(fat) is None or self.sid2geno.get(mot) is None or self.sid2geno.get(kid) is None:
    #                 continue
    #
    #             xt, nt = self.__transmission_count(self.sid2geno[fat], self.sid2geno[mot], self.sid2geno[kid])
    #             Xt = map(operator.add, Xt, xt)
    #
    #         # #(minor allele in affected sib) - #(minor allele in unaffected sib)
    #         for aff, unaff in ped.get_all_DSP():
    #             if self.sid2geno.get(aff) is None or self.sid2geno.get(unaff) is None:
    #                 continue
    #
    #             xs, ns = self.__discordant_count(self.sid2geno[aff], self.sid2geno[unaff])
    #             Xs = map(operator.add, Xs, xs)
    #
    #         dij = [(xt + xs) * 1.0 for (xt, xs) in zip(Xt, Xs)]
    #         Dij.append(dij)
    #         nosib_Dij.append(Xt)
    #     return Dij, nosib_Dij

    def __transmission_count(self, fat, mot, kid):
        """
        transmission count of a trio, for each variant

        return a list of Xt and a list of Nt
        """
        Xt, Nt = [], []
        for (f, m, k) in zip(fat, mot, kid):
            # TODO: ignore the case of mendelian error
            tdt, informative = self.__tdt(f, m, k)
            Xt.append(tdt)
            Nt.append(informative)
        return (Xt, Nt)

    @staticmethod
    def __tdt(f, m, k):
        if f == -9 or m == -9 or k == -9:
            return (0, 0)

        if f == 1 or m == 1:
            # both informative
            if f == 1 and m == 1:
                if k == 0:
                    return (-2, 1)
                elif k == 1:
                    return (0, 1)
                elif k == 2:
                    return (2, 1)

            # only one parent informative
            if f == 1:
                info, ninfo = f, m
            else:
                info, ninfo = m, f

            if k == 0:
                # info pass 0, and ninfo pass 0
                if ninfo == 2:
                    return (0, 0) # mendelian error
                else:
                    return (-1, 1) # info == 0
            elif k == 2:
                # info pass 1, and ninfo pass 1
                if ninfo == 0:
                    return (0, 0) # mendelian error
                else:
                    return (1, 1) # info == 2
            elif k == 1:
                # info pass 1, and ninfo pass 0
                if ninfo == 0:
                    return (1, 1)
                # info pass 1, and ninfo pass 0
                elif ninfo == 2:
                    return (-1, 1)
            else:
                return (0, 0)
        else:
            return (0, 0)

    @staticmethod
    def __discordant_count(aff, unaff):
        """
        discordant count of a DSP, for each variant

        return a list of Xs and a list of Ns
        """
        Xs, Ns = [], []
        for (a, u) in zip(aff, unaff):
            if a == u or a == -9 or u == -9:
                sib, informative = 0, 0
            else:
                sib = a - u
                informative = 1

            Xs.append(sib)
            Ns.append(informative)
        return (Xs, Ns)

    def _rvgdt_statistic_matrix(self):
        """
        RVGDT: Calculate socre matrix Sij (i for pedigree, j for variant)
        # TODO: incorporate covariates (equation 2 in Chen et al 2009)
        """
        has_covariates = False
        if self.pedigree[0].has_covariates is True:
            has_covariates = True

            phenos, covars = [], []
            for ped in self.pedigree:
                for sid, person in ped.sid2subject.iteritems():
                    if sid in self.sid2geno:
                        phenos.append(person.phenotype)
                        covars.append(person.covariates)

            if len(phenos) == 0 or len(phenos) != len(covars):
                has_covariates = False
            else:
                alpha = self.__logit_params(phenos, covars)

        Sij = []
        for ped in self.pedigree:
            Ni = 0 # number of genotyped individuals in pedigree i
            for sid in ped.sid2subject:
                if sid in self.sid2geno:
                    Ni += 1

            if Ni == 0:
                continue

            sij = None
            for aff, unaff in ped.get_discordant_pairs():
                if (aff in self.sid2geno) and (unaff in self.sid2geno):
                    # informative_pedigree = True
                    d = self.__score_diff(self.sid2geno[aff], self.sid2geno[unaff])

                    if has_covariates is True:
                        Zjk = [j - k for j,k in zip(ped.sid2subject[aff].covariates, ped.sid2subject[unaff].covariates)]
                        Zalpha = np.exp(sum(map(operator.mul, Zjk, alpha)))
                        weight = (8.0 / Ni) * Zalpha / pow(1 + Zalpha, 3)
                        d = [i * weight for i in d]

                    sij = d if sij is None else map(operator.add, sij, d)
                    # print d, sij

            if has_covariates is False:
                sij = [s * 1.0 / Ni for s in sij]

            if sij is None:
                continue
            else:
                Sij.append(sij)
        return Sij


    @staticmethod
    def __logit_params(phenos, covars):
        # convert to 0/1
        phenos = [i -1 for i in phenos]
        # add intercept
        covars = [[1] + list(i) for i in covars]
        logit = sm.Logit(phenos, covars)
        try:
            result = logit.fit(disp=0)
        except:
            sys.exit('Unable to infer log odds ratio from phenotype and covariates')
        return result.params[1:]

    @staticmethod
    def __score_diff(aff, unaff):
        """
        genotype score difference between affected and unaffected subject, for each variant
        return a list of integer
        """
        diff = []
        for (a, u) in zip(aff, unaff):
            if a == -9 or u == -9:
                diff.append(0)
            else:
                diff.append(a - u)
        return diff

    def __str__(self):
        mystr = ""
        for ped in self.pedigree:
            mystr += "Pedigree " + ped.fid + "\n"
            for trio in ped.get_all_triad():
                mystr += "Trio {0} + {1} -> {2}\n".format(*trio)
                mystr += "fat: " + " ".join([str(i) for i in self.sid2geno[trio[0]]]) + "\n"
                mystr += "mot: " + " ".join([str(i) for i in self.sid2geno[trio[1]]]) + "\n"
                mystr += "kid: " + " ".join([str(i) for i in self.sid2geno[trio[2]]]) + "\n"
        return mystr


# if __name__ == '__main__':
#     peds, sid2geno = read_pedfile("../example/example.ped")
#
#     weights = [1, 0.5]
#     mafs = [0.005 for i in weights]
#
#     mytest = DisequilibriumTest(peds, sid2geno, mafs, weights)
#     mytest.set_test_options(0.01, False, False, True)
#     print mytest.analyze(adapt_iter=1, max_iter=1)

    # fam2ped = {}
    # sid2geno = {}
    # with open("../example/case.ped") as f:
    #     for line in f.readlines():
    #         items = line.strip().split()
    #         fid = items[0]
    #         info = items[1:6]
    #         if fid not in fam2ped:
    #             fam2ped[fid] = [info]
    #         else:
    #             fam2ped[fid].append(info)
    #
    #         sid = items[1]
    #         geno = [int(i) - 1 for i in items[6:]]
    #         geno = [i+j for i,j in zip(geno[0::2], geno[1::2])]
    #         sid2geno[sid] = geno
    #
    # peds = []
    # for fam in fam2ped:
    #     ped = Pedigree(fam, fam2ped[fam])
    #     peds.append(ped)
    #
    # # for s in sid2geno:
    # #     print s, sid2geno[s]
    # # weights from control
    # weights = None
    # with open("../example/ctrl.ped") as f:
    #     for i, line in enumerate(f.readlines()):
    #         if i % 2 != 0:
    #             continue
    #
    #         items = line.strip().split()
    #         geno = [int(i) - 1 for i in items[6:]]
    #         geno = [i+j for i,j in zip(geno[0::2], geno[1::2])]
    #         if weights is None:
    #             weights = geno
    #         else:
    #             weights = map(operator.add, weights, geno)
    #
    # weights = [[200,2000] for i in weights]
    # mafs = [0.005 for i in weights]
    #
    # mytest = DisequilibriumTest(peds, sid2geno, mafs, weights)
    # print mytest.analyze()