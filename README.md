# tcs-snc2011-data

Dataset of test data for paper:

Akira Terui, GPGCD: An iterative method for calculating approximate GCD of
univariate polynomials.  
Theoretical Computer Science, Volume 479 (Symbolic-Numerical
Algorithms), April 2013, 127-149.  
[doi:10.1016/j.tcs.2012.10.023](http://dx.doi.org/10.1016/j.tcs.2012.10.023)
[arXiv:1207.0630](https://arxiv.org/abs/1207.0630)

## Test 1 (Table 1)

All the programs, scripts and test data are stored under directory
"test1".

### Programs

- Algorithm 1: gpgcd-real-gp.mpl
- Algorithm 2: gpgcd-real.mpl

### Test data

| Test no. | Test polynomial | Test script |
|---|---|---|
| Ex. 1 | testpol-10-10-5-0.1-01.mpl | merom/tcs-snc2011-test-1-1.mw |
| Ex. 2 | testpol-20-20-10-0.1-01.mpl | merom/tcs-snc2011-test-1-2.mw |
| Ex. 3 | testpol-40-40-20-0.1-01.mpl | merom/tcs-snc2011-test-1-3.mw |
| Ex. 4 | testpol-60-60-30-0.1-01.mpl | merom/tcs-snc2011-test-1-4.mw |
| Ex. 5 | testpol-80-80-40-0.1-01.mpl | merom/tcs-snc2011-test-1-5.mw |
| Ex. 6 | testpol-100-100-50-0.1-01.mpl | merom/tcs-snc2011-test-1-6.mw |

## Test 2 (Table 2, Table 3)

All the programs, scripts and test data are stored under directory
"test2".

### Programs

#### Test polynomials generator

- genpol-03.mpl: for generating test polynomials with real
  coefficients
- genpol-complex-02.mpl: for generating test polynomials with complex 
  coefficients

#### Test scripts

##### Test script for each approximate GCD algorithm

- testscript-gpgcd.mpl: for GPGCD
- testscript-stln.mpl: for STLN
- testscript-uvgcd.mpl: for UVGCD

#### GPGCD

- gpgcd-real.mpl: program for polynomials with real coefficients
- gpgcd-complex.mpl: program for polynomials with complex coefficients 

#### STLN

Obtain the following files from Erich Kaltofen's web site into "stln"
directory:

- [con_mulpoly_gcd_2.0.mpl](http://www4.ncsu.edu/~kaltofen/software/manystln/con_mulpoly_gcd_2.0.mpl)
- [subroutine.mpl](http://www4.ncsu.edu/~kaltofen/software/manystln/subroutine.mpl) 
  

#### UVGCD

1. Download ApaTools.zip from Zhonggang Zeng's
   [ApaTools](http://homepages.neiu.edu/~zzeng/apatools.htm) page.
2. Extract ApaTools.zip under "test2" directory.
