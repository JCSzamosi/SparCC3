#!/usr/bin/env python

'''
@author: jonathanfriedman

Module for estimating the correlations in the basis when only compositional data is available.
'''

import warnings
import numpy as np
from numpy import (unravel_index, argmax, ones, corrcoef, cov, r_, 
                   diag, sqrt, where, nan)
from core_methods import to_fractions
from compositional_methods import variation_mat, clr    
from analysis_methods import correlation
try:
    from scipy.stats import nanmedian
except ImportError:
    from numpy import nanmedian
	

def append_indices(excluded,exclude):
    '''
    Append the indx of current excluded value to tuple of previously excluded values.
    '''
    if excluded is None: inds = exclude
    else:                inds = (r_[excluded[0],exclude[0]], r_[excluded[1],exclude[1]])
    return inds
    
def new_excluded_pair(C, th=0.1, previously_excluded=[]):
    '''
    Find component pair with highest correlation among pairs that 
    weren't previously excluded.
    Return the i,j of pair if it's correlaiton >= than th.
    Otherwise return None.
    '''
#    C_temp = abs(C - diag(diag(C)) )
    C_temp = np.triu(abs(C),1) # work only on upper triangle, excluding diagonal
    C_temp[zip(*previously_excluded)] = 0 
    i,j = unravel_index(argmax(C_temp), C_temp.shape) 
    cmax = C_temp[i,j]
    if cmax > th:
        return i,j
    else:  
        return None

def basis_var(f, Var_mat, M, **kwargs):
    '''
    Estimate the variances of the basis of the compositional data x.
    Assumes that the correlations are sparse (mean correlation is small).
    The element of V_mat are refered to as t_ij in the SparCC paper.
    '''
    ## compute basis variances
    try:    M_inv = np.linalg.inv(M)
    except: M_inv = np.linalg.pinv(M)
    V_vec  = Var_mat.sum(axis=1) # elements are t_i's of SparCC paper
    V_base = np.dot(M_inv, V_vec)  # basis variances. 
    ## if any variances are <0 set them to V_min
    V_min  = kwargs.get('V_min', 1e-10)
    V_base[V_base <= 0] = V_min 
    return V_base

def C_from_V(Var_mat, V_base):
    '''
    Given the estimated basis variances and observed fractions variation matrix, 
    compute the basis correlation & covaraince matrices.
    '''
    Vi, Vj = np.meshgrid(V_base, V_base)
    Cov_base = 0.5*(Vi + Vj - Var_mat)
    C_base = Cov_base/ sqrt(Vi) / sqrt(Vj)
    return C_base, Cov_base

def run_sparcc(f, **kwargs):
    '''
    Estimate the correlations of the basis of the compositional data f.
    Assumes that the correlations are sparse (mean correlation is small).
    '''
    th    = kwargs.get('th', 0.1)
    xiter = kwargs.get('xiter', 10)
    ## observed log-ratio variances
    Var_mat = variation_mat(f)
    Var_mat_temp = Var_mat.copy()
    ## Make matrix from eqs. 13 of SparCC paper such that: t_i = M * Basis_Varainces
    D = Var_mat.shape[0] # number of components
    M = ones((D,D)) + diag([D-2]*D) 
    ## get approx. basis variances and from them basis covariances/correlations 
    V_base = basis_var(f, Var_mat_temp, M)
    C_base, Cov_base = C_from_V(Var_mat, V_base)
    ## Refine by excluding strongly correlated pairs
    excluded_pairs = []
    excluded_comp  = np.array([])
    for xi in range(xiter):
        # search for new pair to exclude
        to_exclude = new_excluded_pair(C_base, th, excluded_pairs) #i,j pair, or None
        if to_exclude is None: #terminate if no new pairs to exclude
            break
        # exclude pair
        excluded_pairs.append(to_exclude)
        i,j = to_exclude
        M[i,j] -= 1
        M[j,i] -= 1
        M[i,i] -= 1
        M[j,j] -= 1
        inds = zip(*excluded_pairs)
        Var_mat_temp[inds]   = 0
        Var_mat_temp.T[inds] = 0
        # search for new components to exclude
        nexcluded = np.bincount(np.ravel(excluded_pairs)) #number of excluded pairs for each component
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp      = where(nexcluded>=D-3)[0]
        excluded_comp_new  = set(excluded_comp) - excluded_comp_prev
        if len(excluded_comp_new)>0:
            print excluded_comp
            # check if enough components left 
            if len(excluded_comp) > D-4:
                warnings.warn('Too many component excluded. Returning clr result.')
                return run_clr(f)
            for xcomp in excluded_comp_new:
                Var_mat_temp[xcomp,:] = 0
                Var_mat_temp[:,xcomp] = 0
                M[xcomp,:] = 0
                M[:,xcomp] = 0
                M[xcomp,xcomp] = 1
        # run another sparcc iteration
        V_base = basis_var(f, Var_mat_temp, M)
        C_base, Cov_base = C_from_V(Var_mat, V_base)
        # set excluded components infered values to nans
        for xcomp in excluded_comp:
            V_base[xcomp] = nan
            C_base[xcomp,:] = nan
            C_base[:,xcomp] = nan
            Cov_base[xcomp,:] = nan
            Cov_base[:,xcomp] = nan
    return V_base, C_base, Cov_base

def run_clr(f):
    '''
    Estimate the correlations of the compositional data f.
    Data is transformed using the central log ratio (clr) transform.
    '''
    z        = clr(f)
    Cov_base = cov(z, rowvar=0)
    C_base   = corrcoef(z, rowvar=0)
    V_base   = diag(Cov_base)
    return V_base, C_base, Cov_base
        
def basis_corr(f, method='sparcc', **kwargs):
    '''
    Compute the basis correlations between all components of 
    the compositional data f. 
    
    Parameters
    ----------
    f : array_like
        2D array of relative abundances. 
        Columns are counts, rows are samples. 
    method : str, optional (default 'SparCC')
        The algorithm to use for computing correlation.
        Supported values: SparCC, clr, pearson, spearman, kendall
        Note that the pearson, spearman, kendall methods are not
        altered to account for the fact that the data is compositional,
        and are provided to facilitate comparisons to 
        the clr and sparcc methods.

    Returns
    -------
    V_base: array
        Estimated basis variances.
    C_base: array
        Estimated basis correlation matrix.
    Cov_base: array
        Estimated basis covariance matrix.

    =======   ============ =======   ================================================
    kwarg     Accepts      Default   Desctiption
    =======   ============ =======   ================================================
    th        0<th<1       0.1       exclusion threshold for SparCC.
    xiter     int          10        number of exclusion iterations for SparCC.
    =======   ============ ========= ================================================
    '''    
    method = method.lower()
    k = f.shape[1]
    ## compute basis variances & correlations
    if k<4: 
        raise ValueError, 'Can not detect correlations between compositions of <4 components (%d given)' %k     
    if method == 'clr':
        V_base, C_base, Cov_base = run_clr(f)
    elif method == 'sparcc':
        V_base, C_base, Cov_base = run_sparcc(f, **kwargs)
        tol = 1e-3 # tolerance for correlation range
        if np.max(np.abs(C_base)) > 1 + tol:
            warnings.warn('Sparcity assumption violated. Returning clr result.')
            V_base, C_base, Cov_base = run_clr(f)    
    else:
        raise ValueError, 'Unsupported basis correlation method: "%s"' %method
    return V_base, C_base, Cov_base 

def main(counts, method='SparCC', **kwargs):
    '''
    Compute correlations between all components of counts matrix.
    Run several iterations, in each the fractions are re-estimated, 
    and return the median of all iterations.
    Running several iterations is only helpful with 'dirichlet' 
    normalization method, as with other methods all iterations 
    will give identical results. Thus, if using other normalizations
    set 'iter' parameter to 1.
     
    Parameters
    ----------
    counts : DataFrame
        2D array of counts. Columns are components, rows are samples.
        If using 'dirichlet' or 'pseudo' normalization, 
        counts (positive integers) are required to produce meaningful results, 
        though this is not explicitly checked by the code.  
    method : str, optional (default 'SparCC')
        The algorithm to use for computing correlation.
        Supported values: SparCC, clr, pearson, spearman, kendall
        Note that the pearson, spearman, kendall methods are not
        altered to account for the fact that the data is compositional,
        and are provided to facilitate comparisons to 
        the clr and sparcc methods.

    Returns
    -------
    cor_med: array
        Estimated correlation values.
    cov_med: array
        Estimated covariance matrix if method in {SparCC, clr},
        None otherwise.
              
    =======   ============ =======   ================================================
    kwarg     Accepts      Default   Desctiption
    =======   ============ =======   ================================================
    iter      int          20        number of estimation iteration to average over.
    oprint    bool         True      print iteration progress?
    th        0<th<1       0.1       exclusion threshold for SparCC.
    xiter     int          10        number of exclusion iterations for sparcc.
    norm      str          dirichlet method used to normalize the counts to fractions.
    log       bool         True      log-transform fraction? used if method ~= SparCC/CLR
    =======   ============ ========= ================================================
    '''
    method = method.lower()
    cor_list = []  # list of cor matrices from different random fractions
    var_list = []  # list of cov matrices from different random fractions
    oprint   = kwargs.pop('oprint',True)
    n_iter     = kwargs.pop('iter',20)  # number of iterations 
    norm     = kwargs.pop('norm','dirichlet')
    log      = kwargs.pop('log','True')
    th       = kwargs.setdefault('th',0.1)   # exclusion threshold for iterative sparse algo
    if method in ['sparcc', 'clr']: 
        for i in range(n_iter):
            if oprint: print '\tRunning iteration' + str(i)
            fracs = to_fractions(counts, method=norm)
            v_sparse, cor_sparse, cov_sparse = basis_corr(fracs, method=method, **kwargs)
            var_list.append(np.diag(cov_sparse))
            cor_list.append(cor_sparse)
        cor_array = np.array(cor_list)
        var_med = nanmedian(var_list,axis=0) #median variances
        cor_med = nanmedian(cor_array,axis=0) #median correlations
        x,y     = np.meshgrid(var_med,var_med)
        cov_med = cor_med * x**0.5 * y**0.5
    elif method in ['pearson', 'kendall', 'spearman']:
	n = counts.shape[1]
	cor_array = np.zeros((n_iter, n, n))
        for i in range(n_iter):
            if oprint: print '\tRunning iteration ' + str(i)
            fracs = to_fractions(counts, method=norm)
            if log:
                x = np.log(fracs)
            else:
                x = fracs
            cor_mat, pval = correlation(x, method, axis=0)
            cor_array[i,:,:] = cor_mat
        cor_med = np.median(cor_array, axis=0) #median correlation
        cov_med = None
    return cor_med, cov_med 
        

if __name__ == '__main__':
    ## parse input arguments
    from optparse import OptionParser
    kwargs = {}
    usage  = ('Compute the correlation between components (e.g. OTUs).\n' 
              'By default uses the SparCC algorithm to account for compositional effects.\n' 
              'Correlation and covariance (when applies) matrices are written out as txt files. \n'
              'Counts file needs to be a tab delimited text file where columns are samples and rows are components (e.g. OTUS).\n'
              ' See example/fake_data.txt for an example file.\n' 
              '\n'
              'Usage:   python SparCC.py counts_file [options]\n'
              'Example: python SparCC.py example/fake_data.txt -i 20 --cor_file=example/basis_corr/cor_mat_sparcc.out')
    parser = OptionParser(usage)
    parser.add_option("-c", "--cor_file", dest="cor_file", type = 'str',
                      help="File to which correlation matrix will be written.")
    parser.add_option("-v", "--cov_file", dest="cov_file", type = 'str',
                      help="File to which covariance matrix will be written.")
    parser.add_option("-a", "--algo", dest="algo", default='SparCC',
                      help="Name of algorithm used to compute correlations (SparCC (default) | pearson | spearman | kendall)")
    parser.add_option("-i", "--iter", dest = 'iter', type ='int', default=20,
                      help="Number of inference iterations to average over (20 default).")
    parser.add_option("-x", "--xiter", dest = 'xiter', type ='int', default=10,
                      help="Number of exclusion iterations to remove strongly correlated pairs (10 default).")
    parser.add_option("-t", "--thershold", dest = 'th', type ='float', default=0.1,
                      help= "Correlation strength exclusion threshold (0.1 default).")
    (options, args) = parser.parse_args()
    counts_file     = args[0]
    
    from analysis_methods import basis_corr
    from io_methods import read_txt, write_txt

    kwargs = options.__dict__
    algo     = kwargs.pop('algo')
    cor_file = kwargs.pop('cor_file')
    cov_file = kwargs.pop('cov_file')
    if cor_file is None: cor_file =  'cor_mat_' + algo + '.out' 
    if cov_file is None: cov_file =  'cov_mat_' + algo + '.out' 
        
    print 'reading data'
    counts = read_txt(counts_file)
     
    ## Calculate correlations between components using SparCC
    print 'computing correlations'
    cor, cov = basis_corr(counts, method=algo, **kwargs)
     
    ## write out results
    print 'writing results'
    write_txt(cor, cor_file)
    print 'wrote ' + cor_file
    if cov is not None:
        write_txt(cov, cov_file)
        print 'wrote ' + cov_file
     
    print 'Done!'



