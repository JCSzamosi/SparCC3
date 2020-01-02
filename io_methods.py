'''
Created on Dec 6, 2012

@author: jonathanfriedman
'''
import numpy as np
from pandas.io.parsers import read_table
from Lineages import Lineages
#from pandas.util.decorators import Appender



#@Appender(read_table.__doc__)
def read_txt(file, T=True, lin=None, lin_label='lineage', 
             format='QIIME', verbose=True, **kwargs):
    '''
    Read general delimited file into DataFrame.
    
    This a wrapper around pandas' read_table function which adds
    optional parsing of lineage information, and sets some default
    parameter values.
    
    Note: 
    By default the data is transposed!
    To avoid this behavior set the parameter 'T' to False.
    
    Parameters
    ----------
    file : string 
        Path to input file.  
    T : bool (default True)
        Indicated whether the produced DataFrame will be transposed.
    lin : bool/None (default None)
        Indicated whether lineage information is given in the input file.
        If None, read_txt tries to infer the presence of 
        lineage information automatically
    lin_label : string (default 'lineage')
        Label of the column containing the lineage information.
    format : string (default 'QIIME')
        Format of the lineage information.
        This argument is passed to the Lineage object constructor.
    verbose : bool (default True)
        Indicated whether to print to screen the parsed table stats.
    
    Returns
    -------
    table : DataFrame
        Parsed table.
    lins : Lineages (optional)
        Parsed Lineages object.
        Returned only if lineage information was parsed.
    '''
    kwargs.setdefault('index_col',0)
    temp = read_table(file, **kwargs)
    # try to decide whether lineages are given, if not specified by user
    if lin is None:
        lin = False
        lin_labels = ('lin','lins','lineage','lineages',
                      'taxon','taxa','rdp')
        for c in temp.columns:
            if hasattr(c, 'lower'):
                if c.lower() in lin_labels:
                    lin = True
                    lin_label = c              
    if lin: # parse lins if needed   
        lins = Lineages.from_dict(temp[lin_label], format=format)
        temp = temp.drop(lin_label,axis=1)
    if T: 
        temp = temp.T
        
    s = ['Finished parsing table.',
         'Table dimensions: (%d,%d)' %temp.shape]
    if T:
        s += ['**** Data has been transposed! ****']
    ncol = min(temp.shape[1],3)
    nrow = min(temp.shape[0],3)
    scol = tuple([ncol] + list(temp.columns[:ncol]))
    srow = tuple([nrow] + list(temp.index[:nrow]))
    s +=  [('First %d column labels are :' 
            + ' ,'.join(['%s']*ncol)) %scol,
           ('First %d row labels are :' 
            + ' ,'.join(['%s']*nrow)) %srow]
        
    table = temp
    if lin:
        if verbose: print '\n'.join(s), '\n'
        return table, lins
    else:
        if verbose: print '\n'.join(s), '\n'
        return table

def write_txt(frame, file, T=True, lin=None, lin_label='lineage', **kwargs):
    '''
    Write frame to txt file.
    
    This a wrapper around pandas' to_csv function which adds
    optional writing of lineage information, and sets some default
    parameter values.
    
    Note: 
    By default the data is transposed!
    To avoid this behavior set the parameter 'T' to False.        
    
    Parameters
    ----------
    file : string 
        Path to input file.  
    T : bool (default True)
        Indicated whether the produced DataFrame will be transposed.
    lin : None/None (default None)
        Lineages object to be included in the output file.
    lin_label : string (default 'lineage')
        Label of the column containing the lineage information.    
    '''
    from pandas import Series
    kwargs.setdefault('sep','\t')
    if T: data = frame.T
    else: data = frame
    if lin is not None:
        d = {}
        for i in data.index:
            if i in lin: 
                d[i] = lin[i].lin_str
            else: 
                d[i] = None
        t = Series(d, name=lin_label)
        data = data.join(t)
    data.to_csv(file,**kwargs)

    
if __name__ == '__main__':
    file = 'demo/data/fake_data_lin.counts'
    t, lin = read_txt(file)
#    write_txt(t*1.1, 'temp.txt', lin=lin)
#    print t
#    print read_txt.__doc__
    
