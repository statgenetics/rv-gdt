#!/usr/bin/env python
#
# File: pedigree
#
__author__ = "Zong-Xiao He"
__copyright__ = "Copyright 2015 Zong-Xiao He"
__date__ = "10/30/15"


import itertools
import sys
from random import shuffle


# plink format: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
#      Family ID
#      Individual ID
#      Paternal ID
#      Maternal ID
#      Sex (1=male; 2=female; other=unknown)
#      Phenotype (1=unaffected; 2=affected; 0=unknown)


PHENO_AFFECTED = 2
PHENO_UNAFFECTED = 1
PHENO_MISSING = 0

class Person(object):
    """
    A class for a subject in genetic study.

    Attributes:
        1. sid: subject id (must unique, will be used to get the genotype)
        2. sex:
        3. phenotype: plink pedigree default coding (0 missing / 1 unaffected / 2 affected)
        4. link to parents
        5. is_founder (bool)
        # 6. has_shuffled (bool): a flag indicate whether this subject has been shuffled (to make sure shuffle parent
        #                         before offspring)
    """

    def __init__(self, sid, sex, pheno, covar = None):
        self.sid = sid
        self.covariates = covar
        self.phenotype = pheno
        self.gender = sex
        self.father = None
        self.mother = None
        self.is_founder = False
        # self.has_shuffled = False

    def set_father(self, fat):
        self.father = fat
        return

    def set_mother(self, mot):
        self.mother = mot
        return

    def __str__(self):
        return self.sid


class Pedigree(object):
    """
    A class for a pedigree (all individuals with the same family id)

    Attributes:
        1. fid: family id
        2. sid2subject: subject id to Person object
        3. founders: list of founder subject id
        4. non_founders: list of non founder subject id
    """
    def __init__(self, fid, members):
        self.fid = str(fid)
        self.sid2subject = {}
        self.founders = [] # sid of founders
        self.non_founders = []
        self.affected = []
        self.unaffected = []
        # self.no_phenotype = []
        self.triad = None
        self.dsp = None

        self.has_covariates = False
        self.members = members
        self.__build_pedigree(members)
        return

    def __build_pedigree(self, members):
        # each line: sid, fatid, motid, sex, pheno, [covariates]
        # build nodes
        for m in members:
            if len(m) < 5:
                continue

            sid, sex, pheno = m[0], m[3], m[4]
            if len(m) > 5:
                self.has_covariates = True
                try:
                    covar = [float(i) for i in m[5:]]
                except:
                    sys.exit("Covariates need to be float number!")
            else:
                covar = None

            try:
                sex, pheno = int(sex), int(pheno)
            except Exception as e:
                print(e)
                print "Error: Sex and Phenotype must be integer in pedigree file"
                raise

            if sid not in self.sid2subject:
                self.sid2subject[sid] = Person(sid, sex, pheno, covar)

            if pheno == PHENO_AFFECTED:
                self.affected.append(sid)
            elif pheno == PHENO_UNAFFECTED:
                self.unaffected.append(sid)
            # else:

        # build edges: kid to parent
        for m in members:
            sid, fatid, motid = m[0], m[1], m[2]

            # founder is defined as: both father id and mother id are 0
            if fatid == '0' or motid == '0':
                self.founders.append(sid)
            else:
                fat = self.sid2subject.get(fatid)
                mot = self.sid2subject.get(motid)

                if fat is None or mot is None:
                    # print "Error: father/mother id ({0}/{1}) for individual ({2}) is not found in family ({3})".format(
                    #     fatid, motid, sid, self.fid
                    # )
                    # father mother not found in this family
                    self.founders.append(sid)
                else:
                    self.non_founders.append(sid)
                    self.sid2subject[sid].set_father(fat)
                    self.sid2subject[sid].set_mother(mot)
        return


    def preprocess_rvpdt(self, genotyped):
        """
        RVPDT: break pedigree into nuclear families, discard nuclear family with missing founder
        """
        nfounders = []
        for s in self.founders:
            if s not in genotyped:
                if s in self.affected:
                    self.affected.remove(s)
                if s in self.unaffected:
                    self.unaffected.remove(s)
            nfounders.append(s)
        self.founders = nfounders

        nnon_founders = []
        for s in self.non_founders:
            if s not in genotyped:
                if s in self.affected:
                    self.affected.remove(s)
                if s in self.unaffected:
                    self.unaffected.remove(s)
                continue

            sub = self.sid2subject[s]
            fat, mot = sub.father, sub.mother
            if (fat is None) or (mot is None) or (fat not in genotyped) or (mot not in genotyped):
                if s in self.affected:
                    self.affected.remove(s)
                if s in self.unaffected:
                    self.unaffected.remove(s)
                continue

            nnon_founders.append(s)
        self.nnon_founders = nnon_founders
        return

    def get_triads(self):
        """
        [RVPDT] get all parent - affected kid trios
        """
        if self.triad is not None:
            return self.triad

        self.triad = []
        for s in self.non_founders:
            sub = self.sid2subject[s]
            if sub.phenotype != PHENO_AFFECTED:
                continue

            fat, mot = sub.father, sub.mother
            if fat is not None and mot is not None:
                self.triad.append([fat.sid, mot.sid, sub.sid])
        return self.triad

    def get_DSPs(self):
        """
        [RVPDT] get all discordant sib pairs (same parents, different phenotype)
        """
        if self.dsp is not None:
            return self.dsp

        pat2kids = {}
        self.dsp = []
        # build map: fat_mot -> [kids(with phenotype)]
        for s in self.non_founders:
            sub = self.sid2subject[s]
            if sub.phenotype == PHENO_MISSING:
                continue

            fat, mot = sub.father, sub.mother
            if fat is None or mot is None:
                continue

            pat = '{0}_{1}'.format(fat.sid, mot.sid)
            if pat not in pat2kids:
                pat2kids[pat] = [sub]
            else:
                pat2kids[pat].append(sub)

        # find all discordant sib-pair
        for p in pat2kids:
            kids = pat2kids[p]
            if len(kids) <= 1:
                continue

            for x,y in itertools.combinations(kids, 2):
                if x.phenotype == y.phenotype:
                    continue

                if x.phenotype == PHENO_AFFECTED and y.phenotype == PHENO_UNAFFECTED:
                    self.dsp.append([x.sid, y.sid])
                elif y.phenotype == PHENO_AFFECTED and x.phenotype == PHENO_UNAFFECTED:
                    self.dsp.append([y.sid, x.sid])
                else:
                    continue
        return self.dsp

    def get_discordant_pairs(self):
        """
        [RVGDT] get all discordant pairs within pedigree
        """
        # if self.discordant_pairs is None:
        discordant_pairs = list(itertools.product(self.affected, self.unaffected))

        return discordant_pairs


    def shuffle_phenotype(self):
        sids, phenos = [], []
        for s in self.sid2subject:
            sids.append(s)
            phenos.append(self.sid2subject[s].phenotype)

        shuffle(phenos)
        self.affected, self.unaffected = [], []
        for s,p in zip(sids, phenos):
            self.sid2subject[s].phenotype = p

            if p == PHENO_AFFECTED:
                self.affected.append(s)
            elif p == PHENO_UNAFFECTED:
                self.unaffected.append(s)

        return

    def __str__(self):
        mystr = "Family ID: " + self.fid + "\n"
        mystr += "Founders: \n"
        for f in self.founders:
            mystr += "\t" + f.__str__() + "\n"
        mystr += "Non-founders: \n"
        for f in self.non_founders:
            mystr += "\t" + f.__str__() + \
                     "\t father " + self.sid2subject[f].father.__str__() + \
                     "\t mother " + self.sid2subject[f].mother.__str__() + "\n"
        return mystr


if __name__ == '__main__':
    fam2ped = {}
    with open("../example/multiplex.ped") as f:
        for line in f.readlines():
            items = line.strip().split()
            fid = items[0]
            info = items[1:]
            if fid not in fam2ped:
                fam2ped[fid] = [info]
            else:
                fam2ped[fid].append(info)

    for fam in fam2ped:
        ped = Pedigree(fam, fam2ped[fam])
        print ped
        print "trios: ", ped.get_triads()
        print "DSPs:  ", ped.get_DSPs()
        #print "pairs: ", ped.get_discordant_pairs()

        ped.shuffle_phenotype()
        print "trios: ", ped.get_triads()
        print "DSPs:  ", ped.get_DSPs()