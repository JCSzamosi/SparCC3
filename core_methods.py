'''
Created on Jun 24, 2012

@author: jonathanfriedman 
'''

from pandas import DataFrame as DF
import numpy as np


def _get_axis(given_axis):
    names0 = set([0,'0', 'rows','index','r'])
    names1 = set([1,'1','cols','columns','c'])
    ## to lower case
    if hasattr(given_axis, 'lower'):
        given_axis = given_axis.lower()
    ## get axis
    if given_axis is None:
        return None            
    elif given_axis in names0:
        axis = 0
    elif given_axis in names1:
        axis = 1
    else:
        raise ValueError('Unsupported axis "%s"' %given_axis)  
    return axis

def get_labels(frame,axis):
    '''
    Get the labels of rows/cols.
    '''
    axis = _get_axis(axis)
    return list(frame._get_axis(axis))

def set_labels(frame,labels,axis, inplace=True):
    '''
    Set the labels of rows/cols.
    '''
    axis = _get_axis(axis)
    if isinstance(labels, dict):
        old = frame._get_axis(axis)
        new = [labels[l] for l in old]
    else:
        new = labels
    if inplace:
        frame._set_axis(1-axis,new)
    else:
        fnew = frame.copy()
        fnew._set_axis(1-axis,new)
        return fnew
                   
def _col_numbers(frame, col_lables):
    '''
    Get the numbers of desired cols.
    '''
    ind = [list(frame.columns).index(c) for c in col_lables]
    return ind

def _row_numbers(frame, row_lables):
    '''
    Get the numbers of desired cols.
    '''
    ind = [list(frame.index).index(r) for r in row_lables]
    return ind

def parse_reducer(reducer):
    '''
    Parse the reducer used to combine filters for filter_by_vals.
    '''
    import operator
    reducer = reducer.strip().lower()
    if hasattr(reducer, '__call__'):
        return reducer
    elif isinstance(reducer,str):
        if reducer=='all':
            return operator.and_
        elif reducer=='any':
            return operator.or_
        else:
            raise ValueError, 'Unsupported value for "how" argument: "%s"' %reducer
    else:
        raise TypeError, 'Unsupported "how" argument type'

def filter_by_vals(frame, criteria, axis='cols', verbose=True, how='all',
                   nan_val=None, norm=False):
    '''
    Filter frame, keeping only cols/rows that pass the filtering criteria.
    See pysurvey.util.filters for more information.
        
    Parameters
    ----------
    frame : frame 
        Frame to be filtered
    criteria : filter-parsable/iterable of filter-parsables
        The filtering criteria. 
        Each criterion can be:
            - A triplet of (actor,comperator,value), where the actor extracts the 
              quantity of interest from each row/column and the comperator compares 
              it to the given value and returns a bool.
              Named actors include: 'sum','avg','med', 'var', 'std' and 'presence'.
              A row/col label can also be used to filter by its values. 
              To filter by the values fo a row/col who's label is a named actor, prefix 
              an underscore to it (e.g. '_sum').
              Name comperators include: '==', '!=', '>', '<', '>=', '<=', 'in'.
            - A function that accepts a Series and returns a bool.
            - A pysurvey.util.filters.Filter object.
    axis :  {0 | 1}
        0 : filter rows.
        1 : filter columns.
    verbose : bool (default True)
        Determines whether to print filtering info.
    how : {'all' (default) | 'any' | callable}
        'all'    - Keep row/cols that pass all filtering criteria.
        'any'    - Keep row/cols that pass any of the filtering criteria.
        callable - to be used to reduce the list of bools returned by the filters 
                   for each row/col.
    nan_val : bool/None (default None) 
        Value to be returned by filter if a nan is encountered.
        If None is given, nan are not treated separately.
    norm : bool (default False)
        Indicates whether to normalize the frame before evaluating the filters.
        The filtering itself is always conducted on the unnormalized frame.
        
    Returns
    -------
    filtered: frame
       Filtered frame (new instance). 
    '''
    from pysurvey.util.filters import parse_filters
    axis = _get_axis(axis)
    if norm:
        x = normalize(frame) 
    else:          
        x = frame
    reducer = parse_reducer(how)
    ## create filters
    filters = parse_filters(criteria, nan_val)
    ## find labels to drop   
    selectors = (x.apply(fil, axis=1-axis) for fil in filters)
    selector  = reduce(reducer, selectors)
    drop      = selector[selector==False].index
    ## do filtering
    filtered = frame.drop(drop, axis=axis)
    ## print message
    if verbose:
        axis_s = {0:'rows',1:'columns'}
        s = ['Dropped %d %s' %(len(drop),axis_s[axis]),
             'Resulting size is (%d,%d)' %filtered.shape]
        print '\n'.join(s) +'\n'
    return filtered

def keep(frame, n, criterion='sum', axis=0, which='first', sort=True):
    '''
    Create a new frame with only the n most extreme rows/cols.
    
    -------- NO UNITTEST ---------
    
    Parameters
    ----------
    frame : frame 
        Frame to be filtered
    n : int
        Number of row/cols to be kept.
    criterion : {'sum' (default) | 'avg' | 'med' | 'std' | 'presence' | 'var' | label | callable}
        Criterion by which the row/columns will be ordered.
        See pysurvey.util.filters.parse_actor for more information.
    axis :  {0 | 1}
        0 : keep only n rows.
        1 : keep only n cols.
    which : {'first' (default) | last}
        Indicates whether to keep the first or last n elements after sorting by criterion.
    sort : bool (default False)
        Indicates whether to sort the kept n rows/cols by the given criterion,
        or retain the order in which they appear in the given frame.
        
    Returns
    -------
    filtered: frame
       Filtered frame (new instance). 
    '''
    from pysurvey.util.filters import parse_actor
    axis = _get_axis(axis)
    if   axis == 1: data = frame
    elif axis == 0: data = frame.T
    f = parse_actor(criterion)
    
#    biggest = kwargs.get('biggest', True) # if true return the n cols with the biggest values for criterion, else return the n cols with the smallest values.
    temp = data.apply(f)
    temp.sort()
    temp = temp[::-1]
    which = which.strip().lower()
    if which == 'first':
        inds = temp.index[:n]
    elif which == 'last': 
        inds = temp.index[-n:]
    else:
        raise ValueError, "Unsupported value for 'which' parameter: %s" %which      
    filtered    = data.filter(items=inds)
    if not sort: filtered = filtered.reindex_like(data).dropna(how='all', axis=axis)
    if axis == 0: filtered = filtered.T
    return filtered
            
def vals_by_keys(frame, key_pairs):
    '''
    Return a list of values corresponding to key_pairs.
    Inputs:
        key_pairs = [list] each element = [col_key, row_key].
    Outputs:
        vals = [list] values for each pair in key_pairs, in corresponding order.
    '''
    vals = map(lambda pair: frame[pair[0]][pair[1]], key_pairs)
    return vals

def to_binary(frame, th=0):
    '''
    Discretize matrix s.t. matrix[matrix > th] = 1, matrix[matrix <= th] = 0. 
    Return new instance.
    '''
    bin = frame.copy()
    ind = frame > th
    bin[ind]  = 1
    bin[-ind] = 0
    return bin


#-------------------------------------------------------------------------------
# Methods for counts data

def normalize(frame, axis=0):
    '''
    Normalize counts by sample total.
    
    Parameters
    ----------
    axis : {0, 1}
        0 : normalize each row
        1 : normalize each column

    Returns new instance of same class as input frame.
    '''
    axis = _get_axis(axis)
    if isinstance(frame, DF):
        return frame.apply(lambda x:1.*x/x.sum(), axis=1-axis)
    else:
        return np.apply_along_axis(lambda x:1.*x/x.sum(), 1-axis, frame)
    
def to_fractions(frame, method='dirichlet', p_counts=1, axis=0):
    '''
    Covert counts to fraction using given method.
    
    Parameters
    ----------
    method : string {'dirichlet' (default) | 'normalize' | 'pseudo'}
        dirichlet - randomly draw from the corresponding posterior 
                    Dirichlet distribution with a uniform prior.
                    That is, for a vector of counts C, 
                    draw the fractions from Dirichlet(C+1). 
        normalize - simply divide each row by its sum.
        pseudo    - add given pseudo count (defualt 1) to each count and
                    do simple normalization.
    p_counts : int/float (default 1)
        The value of the pseudo counts to add to all counts.
        Used only if method is dirichlet
    axis : {0 | 1}
        0 : normalize each row.
        1 : normalize each column.
    
    Returns
    -------
    fracs: frame/array
        Estimated component fractions.
        Returns new instance of same class as input frame.
    '''
    axis = _get_axis(axis)
    if method is 'normalize': 
        fracs = normalize(frame, axis)
        return fracs
    
    ## if method is not normalize, get the pseudo counts (dirichlet prior)
    from numbers import Number
    if not isinstance(p_counts, Number):
        p_counts = np.asarray(p_counts)
        
    if method is 'pseudo': 
        fracs = normalize(frame+p_counts, axis)
    elif method is 'dirichlet':
        from numpy.random.mtrand import dirichlet
        def dir_fun(x):
            a = x+p_counts
            f = dirichlet(a)
            return f
        if isinstance(frame, DF):
            fracs = frame.apply(dir_fun, 1-axis)
        else:
            fracs = np.apply_along_axis(dir_fun, 1-axis, frame)
    else:
        raise ValueError, 'Unsupported method "%s"' %method
    return fracs

def rarefy(frame,n, replace=False, remove_shallow=None):
    '''
    Down-sample all rows to have exactly n counts in total for each row.
    if remove_shallow, samples with less than n total counts are excluded.
    
    Parameters
    ----------
    n : int
        Rows will be down-sampled to this total number of counts.
    replace : bool (default False)
        Indicates whether sampling is done with or without replacement.
    remove_shallow : bool/None (default None)
        Indicates whether to remove rows that have less than n total counts to 
        begin with. 
        If None is given, remove_shallow is set to be False for sampling with replacement
        and True for sampling without replacement.
        If remove_shallow is set to False, and sampling is without replacement, 
        rows that have less than the desired total-number of counts are left unchanged.
          
    Returns
    -------
    deep_rarefied: frame
        Rarefied frame (new instance).
    '''
    ## decide whether to remove 'shallow' samples
    if remove_shallow is None: 
        remove_shallow = not replace
    if remove_shallow:
        deep = filter_by_vals(frame, ('sum','>=', n), axis='rows')
    else: 
        deep = frame
    deep_rarefied = deep.copy()
    ## perform rarefaction
    if replace:
        from numpy.random.mtrand import multinomial
        def draw(x):
            p = x/float(x.sum())
            f = 1.*multinomial(n,p)
            return f
    else:
        from numpy.random import rand
        def draw(x):
            k = len(x)
            nt = x.sum()
            if nt < n:
                return x
            new = np.zeros(k)
            counts = 1.*x
            for j in xrange(n):
                p = counts/nt
                i = np.where((p.cumsum() - rand())>0)[0][0]
                nt-=1
                counts[i]-=1
                new[i]+=1
            return new
    deep_rarefied = (deep_rarefied.T.apply(draw)).T
    return deep_rarefied

def group_taxa(frame, lins, level='p', best=True):
    '''
    Return a new instance with cols corresponding to counts aggregated at the 
    desired taxonomic level (e.g. phylum).
    OTUs that are missing from lin are not accounted for.
    OTUs that are not assigned at desired level are aggregated into the 'unassigned' row.

    Parameters
    ----------
    lins : Lineages
        Lineage info of OTUs in frame.
    level : str {'k' | 'p' (default) | 'c' | 'o' | 'f' | 'g' | 's'}
        Desired taxonomic level of aggregation
    best : bool (default True)
        Indicates whether to return the best assigned taxonomy 
        (at the desired level or above), or return the taxonomy at the desired level,
        even if it is unassigned. 
          
    Returns
    -------
    Grouped frame (new instance).
    '''
    old    = frame               
    new    = old.filter(items = []) # create new object of same class as frame, with same samples but now otus.
    taxa   = set(lins.get_assignments(level, best=best)) # set of all taxa present in lin.
    for t in taxa:
        otus = lins.get_ids(level, t, best=best)
        temp = old.filter(items = otus) # matrix with only otus of given taxa
        new[t] = temp.sum(axis = 1)
    return new.dropna(axis=0)
        

if __name__ == '__main__':
    rows = ['r1', 'r0', 'r2', 'r3']
    cols = ['c0', 'c1', 'c2']
    metac = DF([[np.nan,'big'],
            ['Entero','small'],
            ['Blautia','tiny']], 
           columns=['name', 'Size'],
           index=cols)
    mat = np.array([[2., np.NAN,1], 
                    [1, 3, 2],
                    [10, 15,3],
                    [0,0,1]])
    df = DF(mat, index=rows, columns=cols)
#    print df,'\n'
#    print filter_by_vals(df,[('sum','<=',3),('presence','>',1)],axis='rows'),'\n'   
    print metac, '\n'
    actor = lambda x: x['Size']
    filter1 = lambda x: isinstance(x['name'], str)
    filter2 = (actor,'in',['big','tiny'])
    filter3 = ('Size','in',['big','tiny'])
    filter4 = ('name','in',['Entero','Blautia'])
    print filter_by_vals(metac, filter1, axis=0),'\n'
    print filter_by_vals(metac, filter2, axis=0),'\n'
    print filter_by_vals(metac, filter3, axis=0),'\n'
    print filter_by_vals(metac,[filter1,filter2], axis=0),'\n'
    print filter_by_vals(metac, filter4, axis=0, nan_val=True),'\n'
    
    df = DF([[1,3,2],[4,6,5]], columns=['a','b','c'], index=['r1','r2'])
    print df,'\n'
    print rarefy(df,7, replace=True)


