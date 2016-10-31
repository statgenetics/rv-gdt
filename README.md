# RV-GDT manual

## Installation

### Requirements

+ python: 2.7


+ numpy: >= 1.11.0

+ spicy: >= 0.16.1

  or

+ Anaconda: >= 2.3

The RV-GDT package is free and available on github.  Run the following commands to download and install the RV-GDT.

``` shell
git checkout https://github.com/hezx/RV-GDT.git
cd RV-GDT
python setup.py install 
```

If the program is installed correctly, you will see program options using the following command:

```shell
rvgdt --help
```



## Input Format

### Genotype File

The genotype file gives the genotype information of each subject per line. No header is needed in the genotype file. The first column is the subject id (which should be the same as the subject id in pedigree files), and the following columns is the number of minor allele on each variant sites (0/1/2 coding and -9 is missing), which is separated by a space or tab. An example of the genotype file is given below

```
11000.fa -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11000.mo -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11000.p1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11000.s1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11001.fa -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11001.mo -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11001.p1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11002.fa -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11002.mo -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
11002.p1 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 -9
```



### Pedigree File

The pedigree file is a white-space (space or tab) delimited file without header. The first six columns are mandatory:

+ Family ID     
+ Individual ID: must be unique within the family     
+ Paternal ID: 0 if not available 
+ Maternal ID: 0 if not available 
+ Sex:  1=male, 2=female   
+ Phenotype: 1=unaffected, 2=affected

Any data given from column 7 will be considered as covariates. An example pedigree file is given below:

```
11000 11000.fa 0 0 1 1
11000 11000.mo 0 0 2 1
11000 11000.p1 11000.fa 11000.mo 1 2
11000 11000.s1 11000.fa 11000.mo 2 1
11001 11001.fa 0 0 1 1
11001 11001.mo 0 0 2 1
11001 11001.p1 11001.fa 11001.mo 1 2
11002 11002.fa 0 0 1 1
11002 11002.mo 0 0 2 1
11002 11002.p1 11002.fa 11002.mo 2 2
```

### Optional Files

#### Ungenotyped Subjects

The full pedigree structure needs to be given in the pedigree file, even if some subjects don't have genotype information. The subjects who are not genotyped are given in a separate single-column file (no header), in which each line is the subject id.

#### Weight File

In case you want to weight each variant differently. The weights can be given in a single-column file (no header), in which each line is the weight for the corresponding variant site (the order should be the same as the order of variant sites in genotype file).

## Options

```
positional arguments:
  proj                  Name of the project

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --geno GENO           genotype file
  --ped PED             pedigree file
  --weight WEIGHT       weights for variants
  --maf_cutoff MAF_CUTOFF
                        MAF cutoff, MAF was calculated in founders
  --adapt_iter ADAPT_ITER
                        Adaptive permutation
  --max_iter MAX_ITER   Max number of iterations
  --alpha ALPHA         Alpha level
  --unweighted          include unweighted test in analysis
  --use_vt              using vt approach in analysis
  --remove_subject REMOVE_SUBJECT
                        remove given subjects
```

Example commands are shown below:

```shell
rvgdt test --geno ./example/rvgdt_test.geno --ped ./example/rvgdt_test.ped --max_iter 100

rvgdt test --geno ./example/rvgdt_test.geno --ped ./example/rvgdt_test_covariates.ped --max_iter 100 --weight ./example/rvgdt_test.weights
```

The output is given in the ${proj}.rvgdt_output file. 

# Questions

If you have any further questions, please fell free to send us an email:

+ Zongxiao He: zongxiah@bcm.edu
+ Suzanne M. Leal: steal@bcm.edu

