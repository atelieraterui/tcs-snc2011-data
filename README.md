# Dataset for the GPGCD Algorithm

Dataset of test data for paper:

Akira Terui, GPGCD: An iterative method for calculating approximate GCD of
univariate polynomials.  
Theoretical Computer Science, Volume 479 (Symbolic-Numerical
Algorithms), April 2013, 127-149.  
[doi:10.1016/j.tcs.2012.10.023](http://dx.doi.org/10.1016/j.tcs.2012.10.023)
[arXiv:1207.0630](https://arxiv.org/abs/1207.0630)

## Test 1 (Table 1)

All files are stored under directory "test1".

### Test data

| Test no. | Test polynomial | Test script and log |
|---|---|---|
| Ex. 1 | testpol-10-10-5-0.1-01.mpl | merom/tcs-snc2011-test-1-1.mw |
| Ex. 2 | testpol-20-20-10-0.1-01.mpl | merom/tcs-snc2011-test-1-2.mw |
| Ex. 3 | testpol-40-40-20-0.1-01.mpl | merom/tcs-snc2011-test-1-3.mw |
| Ex. 4 | testpol-60-60-30-0.1-01.mpl | merom/tcs-snc2011-test-1-4.mw |
| Ex. 5 | testpol-80-80-40-0.1-01.mpl | merom/tcs-snc2011-test-1-5.mw |
| Ex. 6 | testpol-100-100-50-0.1-01.mpl | merom/tcs-snc2011-test-1-6.mw |

### Program files

- gpgcd-real-gp.mpl: Algorithm 1
- gpgcd-real.mpl: Algorithm 2 


## Test 2 (Table 2, Table 3)

All files are stored under directory "test2".

### Test data for Table 2

| Test no. | Test polynomial | Test polynomial generator |  Test script and log | Saved test results (GPGCD\|STLN\|UVGCD) |
| --- | --- | --- | --- | --- |
| Ex. 1 | testpol-10-10-5-0.1-02.mpl | gentestpol-10-10-5-0.1-02.mw | merom-real/test-real-10-10-5-0.1-02.mw | merom-real/testdata-10-10-5-0.1-02-(gpgcd\|uvgcd).mpl, testdata-10-10-5-0.1-02-stln-01.mpl (for STLN) |
| Ex. 2 | testpol-20-20-10-0.1-02.mpl | gentestpol-20-20-10-0.1-02.mw | merom-real/test-real-20-20-10-0.1-02.mw | merom-real/testdata-20-20-10-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 3 | testpol-30-30-15-0.1-02.mpl | gentestpol-30-30-15-0.1-02.mw | merom-real/test-real-30-30-15-0.1-02.mw | merom-real/testdata-30-30-15-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 4 | testpol-40-40-20-0.1-02.mpl | gentestpol-40-40-20-0.1-02.mw | merom-real/test-real-40-40-20-0.1-02.mw | merom-real/testdata-40-40-20-0.1-02-(gpgcd\|stln\|uvgcd).mpl | 
| Ex. 5 | testpol-50-50-25-0.1-02.mpl | gentestpol-50-50-25-0.1-02.mw | merom-real/test-real-50-50-25-0.1-02.mw| merom-real/testdata-50-50-25-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 6 | testpol-60-60-30-0.1-02.mpl | gentestpol-60-60-30-0.1-02.mw | merom-real/test-real-60-60-30-0.1-02.mw | merom-real/testdata-60-60-30-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 7 | testpol-70-70-35-0.1-02.mpl | gentestpol-70-70-35-0.1-02.mw | merom-real/test-real-70-70-35-0.1-02.mw | merom-real/testdata-70-70-35-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 8 | testpol-80-80-40-0.1-02.mpl | gentestpol-80-80-40-0.1-02.mw | merom-real/test-real-80-80-40-0.1-02.mw | merom-real/testdata-80-80-40-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 9 | testpol-90-90-45-0.1-02.mpl | gentestpol-90-90-45-0.1-02.mw | merom-real/test-real-90-90-45-0.1-02.mw | merom-real/testdata-90-90-45-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 10 | testpol-100-100-50-0.1-01.mpl | gentestpol-100-100-50-0.1-02.mw | merom-real/test-real-100-100-50-0.1-02.mw | merom-real/testdata-100-100-50-0.1-02-(gpgcd\|stln\|uvgcd).mpl |

### Equation number in the test log for data in Table 2

Each data in Table 2 appears in test log file ("Test script and log" in the above) with corresponding equation number.

| Perturbation ||| Time ||| #iterations ||
| --- | --- | --- | --- | --- | --- | --- | --- |
| STLN | UVGCD | GPGCD | STLN | UVGCD | GPGCD | STLN | GPGCD |
| (3.2.5) | (4.2.4) | (2.2.11)  | (3.2.2) | (4.2.2) | (2.2.9) | (3.2.4) | (2.2.10) |

### Test data for Table 3

| Test no. | Test polynomial | Test polynomial generator |  Test script and log | Saved test results (GPGCD\|STLN\|UVGCD) |
|---|---|---|---|---|
| Ex. 1 | testpol-complex-10-10-5-0.1-02.mpl | gentestpol-complex-10-10-5-0.1-02.mw | merom-complex/test-complex-10-10-5-0.1-02.mw | merom-complex/testdata-complex-10-10-5-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 2 | testpol-complex-20-20-10-0.1-02.mpl | gentestpol-complex-20-20-10-0.1-02.mw | merom-complex/test-complex-20-20-10-0.1-02.mw | merom-complex/testdata-complex-20-20-10-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 3 | testpol-complex-30-30-15-0.1-02.mpl | gentestpol-complex-30-30-15-0.1-02.mw | merom-complex/test-complex-30-30-15-0.1-02.mw | merom-complex/testdata-complex-30-30-15-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 4 | testpol-complex-40-40-20-0.1-02.mpl | gentestpol-complex-40-40-20-0.1-02.mw | merom-complex/test-complex-40-40-20-0.1-02.mw | merom-complex/testdata-complex-40-40-20-0.1-02-(gpgcd\|stln\|uvgcd).mpl | 
| Ex. 5 | testpol-complex-50-50-25-0.1-02.mpl | gentestpol-complex-50-50-25-0.1-02.mw | merom-complex/test-complex-50-50-25-0.1-02.mw| merom-complex/testdata-complex-50-50-25-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 6 | testpol-complex-60-60-30-0.1-02.mpl | gentestpol-complex-60-60-30-0.1-02.mw | merom-complex/test-complex-60-60-30-0.1-02.mw | merom-complex/testdata-complex-60-60-30-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 7 | testpol-complex-70-70-35-0.1-02.mpl | gentestpol-complex-70-70-35-0.1-02.mw | merom-complex/test-complex-70-70-35-0.1-02.mw | merom-complex/testdata-complex-70-70-35-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 8 | testpol-complex-80-80-40-0.1-02.mpl | gentestpol-complex-80-80-40-0.1-02.mw | merom-complex/test-complex-80-80-40-0.1-02.mw | merom-complex/testdata-complex-80-80-40-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 9 | testpol-complex-90-90-45-0.1-02.mpl | gentestpol-complex-90-90-45-0.1-02.mw | merom-complex/test-complex-90-90-45-0.1-02.mw | merom-complex/testdata-complex-90-90-45-0.1-02-(gpgcd\|stln\|uvgcd).mpl |
| Ex. 10 | testpol-complex-100-100-50-0.1-01.mpl | gentestpol-complex-100-100-50-0.1-02.mw | merom-complex/test-complex-100-100-50-0.1-02.mw | merom-complex/testdata-complex-100-100-50-0.1-02-(gpgcd\|stln\|uvgcd).mpl |

### Equation number in the test log for data in Table 3

Each data in Table 3 appears in test log file ("Test script and log" in the above) with corresponding equation number.

| Perturbation ||| Time ||| #iterations ||
| --- | --- | --- | --- | --- | --- | --- | --- |
| STLN | UVGCD | GPGCD | STLN | UVGCD | GPGCD | STLN | GPGCD |
| (3.2.5) | (4.2.4) | (2.2.5)  | (3.2.2) | (4.2.2) | (2.2.2) | (3.2.4) | (2.2.4) |

### Program files

#### GPGCD (by the author)

Program files for the GPGCD algorithm.

- gpgcd-real.mpl: program for polynomials with real coefficients.
- gpgcd-complex.mpl: program for polynomials with complex coefficients.

#### STLN (by Erich Kaltofen, et al.)

Program files for the STLN-based approximate GCD algorithm.

- con_mulpoly_gcd_2.0-01.mpl: program for polynomials with real
  coefficients.
- con_mulpoly_gcd_2.0-02.mpl: program for polynomials with complex
  coefficients.
  
To obtain the above files, follow the instructions below:

1. Download the following files from Kaltofen's web site:
    - [con_mulpoly_gcd_2.0.mpl](http://www4.ncsu.edu/~kaltofen/software/manystln/con_mulpoly_gcd_2.0.mpl)
    - [subroutine.mpl](http://www4.ncsu.edu/~kaltofen/software/manystln/subroutine.mpl) 
1. Place 2 files in "stln" directory.
1. Apply patches in "stln" directory as follows:

```
cd stln/
patch -c con_mulpoly_gcd_2.0.mpl con_mulpoly_gcd_2.0-01.mpl.patch -o con_mulpoly_gcd_2.0-01.mpl
patch -c con_mulpoly_gcd_2.0.mpl con_mulpoly_gcd_2.0-02.mpl.patch -o con_mulpoly_gcd_2.0-02.mpl
```
Patches have been written by the author.

#### UVGCD (by Zhonggang Zeng)

Program files for the UVGCD algorithm, which can be obtained by the
following instructions:

1. Download ApaTools.zip from Zen's
   [ApaTools](http://homepages.neiu.edu/~zzeng/apatools.htm) page.
2. Extract ApaTools.zip under "test2" directory (so that a directory
   "ApaTools" appears right under "test2" directory).

#### Test polynomial generators

Scripts used for generating test polynomials.

- genpol-03.mpl: for generating test polynomials with real
  coefficients
- genpol-complex-02.mpl: for generating test polynomials with complex 
  coefficients

#### Test scripts

Scripts driven for tests of each approximate GCD algorithm.

- testscript-gpgcd.mpl: for GPGCD.
- testscript-stln.mpl: for STLN.
- testscript-uvgcd.mpl: for UVGCD.

#### Utilities

Utilities used in test scripts.

- largeElementIndex.mpl: for array of numbers, find numbers in the
array which are larger than given criterion and return its index.
- smallElementIndex.mpl: for array of numbers, find numbers in the
array which are smaller than given criterion and return its index.

-----
## Tables 4--9

Under preparation.
