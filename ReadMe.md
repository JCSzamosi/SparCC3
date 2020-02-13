SparCC for Python3
==================

This is a port of SparCC by Jonathan Friedman ([paper
here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687),
original [bitbucket repo
here](https://bitbucket.org/yonatanf/sparcc/src/default/)) to Python3, including
a few fixes. So far I have fixed it to the point that the example data included
with the scripts runs without raising any errors. I would be thrilled to accept
pull requests that do more maintenance. Importantly, only the SparCC correlation 
method is working at the moment. None of the other methods (Pearson, Spearman, etc.)
are correctly implemented.

Original Documentation
----------------------

*The below documentation is adapted from the [original SparCC bitbucket
repository](https://bitbucket.org/yonatanf/sparcc/src/default/).*

`SparCC` is a python module for computing correlations in compositional data
(16S, metagenomics, etc').

Detailed information about the algorithm can be found in the accompanying
[publication](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687).  

Questions, comments, complaints and praise should be send to yonatanf@mit.edu



### Usage Notes:

* Scripts in the root SparCC directory can be called from the terminal
  command-line either by explicitly calling python (as is done in the usage
  examples below), or simply as an executable. The latter will require having
  execution permission for these file (e.g. chmod +x SparCC.py).

* Help for any one for the scripts in the root SparCC directory is available by
  typing 'python [script_name] - h' in the command line. e.g.: 

```
   python SparCC.py -h .
```

* SparCC is implemented in pure python and requires a working version of python
  (=>3.7, tested with 3.7.3), numpy (tested with version 1.17.4), and pandas
  (tested with version 0.25.3).


### Usage example:

* The following lists the commands required for analyzing the included 'fake'
  dataset using the SparCC package, and generating all the files present in the
  subfolders of the example folder.

* The fake dataset contains simulated abundances of 50 otus in 200 samples,
  drawn at random from a multinomial log-normal distribution. The true basis
  correlations used to generate the data are listed in 'true_basis_cor.txt' in
  the example folder.

* Note that otu 0 is very dominant, and thus, using Pearson or Spearman
  correlations, appears to be negatively correlated with most other OTUs, though
  it is in fact not negatively correlated with any OTU.


### Correlation Calculation:

First, we'll quantify the correlation between all OTUs, using SparCC, Pearson,
and Spearman correlations:

```
python SparCC.py example/fake_data.txt -i 5 --cor_file=example/basis_corr/cor_sparcc.out
python SparCC.py example/fake_data.txt -i 5 --cor_file=example/basis_corr/cor_pearson.out -a pearson
python SparCC.py example/fake_data.txt -i 5 --cor_file=example/basis_corr/cor_spearman.out -a spearman
```

### Pseudo p-value Calculation:

Calculating pseudo p-values is done via a bootstrap procedure.
First make shuffled (w. replacement) datasets:

```
python MakeBootstraps.py example/fake_data.txt -n 5 -t permutation_#.txt -p example/pvals/
```

This will generate 5 shuffled datasets, which is clearly not enough to get
meaningful p-values, and is used here for convenience.  A more appropriate
number of shuffles should be at least a 100, which is the default value. 

Next, you'll have to run SparCC on each of the shuffled data sets.  Make sure to
use the exact same parameters which you used when running SparCC on the real
data, name all the output files consistently, numbered sequentially, and with a
'.txt' extension.

```
python SparCC.py example/pvals/permutation_0.txt -i 5 --cor_file=example/pvals/perm_cor_0.txt
python SparCC.py example/pvals/permutation_1.txt -i 5 --cor_file=example/pvals/perm_cor_1.txt
python SparCC.py example/pvals/permutation_2.txt -i 5 --cor_file=example/pvals/perm_cor_2.txt
python SparCC.py example/pvals/permutation_3.txt -i 5 --cor_file=example/pvals/perm_cor_3.txt
python SparCC.py example/pvals/permutation_4.txt -i 5 --cor_file=example/pvals/perm_cor_4.txt
```

Above I'm simply called SparCC 5 separate times. However, it is much more
efficient and convenient to write a small script that automates this, and
submits these runs as separate jobs to a cluster (if one is available to you.
Otherwise, this may take a while to run on a local machine...).

Now that we have all the correlations computed from the shuffled datasets, we're
ready to get the pseudo p-values.  Remember to make sure all the correlation
files are in the same folder, are numbered sequentially, and have a '.txt'
extension.  The following will compute both one and two sided p-values.

```
 python PseudoPvals.py example/basis_corr/cor_sparcc.out example/pvals/perm_cor_#.txt 5 -o example/pvals/pvals.one_sided.txt -t one_sided
 python PseudoPvals.py example/basis_corr/cor_sparcc.out example/pvals/perm_cor_#.txt 5 -o example/pvals/pvals.one_sided.txt -t two_sided
 ```



LICENSE
===================

The MIT License (MIT)

Copyright (c) 2018-2020 Jonathan Friedman and Eric Alm

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



