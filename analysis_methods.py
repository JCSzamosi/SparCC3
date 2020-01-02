'''
Created on Jun 24, 2012

@author: jonathanfriedman 
'''

from pandas import DataFrame as DF
from core_methods import _get_axis
import numpy as np

def basis_corr(frame, algo='SparCC', **kwargs):
    '''
    Compute correlations between all columns of a counts frame.
    This is a wrapper around pysurvey.analysis.basis_correlations.main
        
    Parameters
    ----------
    counts : array_like
        2D array of counts. Columns are components, rows are samples. 
    method : str {SparCC (default)| clr| pearson| spearman| kendall}
        The algorithm to use for computing correlation.

    Returns
    -------
    cor_med: frame
        Estimated correlation matrix.
        Labels are column labels of input frame.
    cov_med: frame/None
        If method in {SparCC, clr} : Estimated covariance matrix.
        Labels are column labels of input frame. 
        Otherwise: None.
              
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
    import SparCC
    comps  = frame.columns
    cor_med, cov_med = SparCC.main(frame, algo=algo, **kwargs)
    print cor_med.shape
    cor = DF(cor_med, index=comps, columns=comps)
    if cov_med is None:
        cov = None
    else:
        cov  = DF(cov_med, index=comps, columns=comps)
    return cor, cov  

def correlation(frame, method='pearson', axis=0):
    '''
    Calculate the correlation between all rows/cols.
    Return frames of correlation values and p-values.
    
    Parameters
    ----------
    frame : DataFrame
        Frame containing data.
    method : {pearson (default) | spearman | kendall}
        Type of correlations to be computed
    axis : {0, 1}
        - 0 - Compute correlation between columns
        - 1 - Compute correlation between rows
    
    Returns
    -------
    c : frame
        DataFrame of symmetric pairwise correlation coefficients.
        Labels are the rows/column labels of the input frame.
    p : frame
        DataFrame of p-values associated with correlation values.
        Labels are the rows/column labels of the input frame.
    ''' 
    import scipy.stats as stats
    axis = _get_axis(axis)
    method = method.lower()
    if method not in set(['pearson', 'kendall', 'spearman']): 
        raise ValueError('Correlation of method %s is not supported.' %method)
    if method == 'spearman' : 
        c_mat, p_mat = stats.spearmanr(frame.values, axis=axis)
        if not np.shape(c_mat):
            c_mat = np.array([[1, c_mat],[c_mat,1]])
            p_mat = np.array([[1, p_mat],[p_mat,1]])
        labels = frame._get_axis(1-axis)
        c = DF(c_mat, index=labels, columns=labels)
        p = DF(p_mat, index=labels, columns=labels)
    else:
        if   method == 'pearson': corr_fun = stats.pearsonr
        elif method == 'kendall': corr_fun = stats.kendalltau
        if   axis == 0: data = frame.T
        elif axis == 1: data = frame
        mat = data.values
        row_labels = data.index
        n = len(row_labels)
        c_mat = np.zeros((n, n))
        p_mat = np.zeros((n, n))
        for i in xrange(n):
            for j in xrange(i, n):
                if i == j: 
                    c_mat[i][i] = 1
                    p_mat[i][i] = 1
                    continue
                c_temp, p_temp = corr_fun(mat[i, :], mat[j, :])
                c_mat[i][j] = c_temp
                c_mat[j][i] = c_temp
                p_mat[i][j] = p_temp
                p_mat[j][i] = p_temp
        c = DF(c_mat, index=row_labels, columns=row_labels)
        p = DF(p_mat, index=row_labels, columns=row_labels)
    return c, p


#-------------------------------------------------------------------------------
# Misc.                    
def permute_w_replacement(frame, axis=0):
    '''
    Permute the frame values across the given axis.
    Create simulated dataset were the counts of each component (column)
    in each sample (row), are randomly sampled from the all the 
    counts of that component in all samples.
    
    Parameters
    ----------
    frame : DataFrame
        Frame to permute.
    axis : {0, 1}
        - 0 - Permute row values across columns
        - 1 - Permute column values across rows    
    
    Returns
    -------
    Permuted DataFrame (new instance).
    '''
    from numpy.random import randint 
    axis = 1-_get_axis(axis)
    s = frame.shape[axis]
    fun = lambda x: x.values[randint(0,s,(1,s))][0]
    perm = frame.apply(fun, axis=axis)
    return perm




if __name__ == '__main__':
    pass
    
    
