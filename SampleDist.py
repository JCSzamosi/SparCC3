#!/usr/bin/env python

'''
Created on Jun 20, 2011

@author: jonathanfriedman

Requires the scipy.cluster.hierarchy module!
'''


from lib.SurveyMatrix import Survey_matrix as SM

    
def kwargs_callback(option, opt, value, parser,**kwargs):
    d = kwargs['d']
    d[option.dest] = value
    return d


def Run(counts_file, metric = 'JSsqrt', **kwargs):
    '''
    Compute the pairwise distance matrix between all sites and write it out as txt file.
    '''
    ## read counts data
    temp   = SM()
    counts = temp.from_file(counts_file)
    ## compute sample distances
    fracs = counts.to_fractions('normalize')
    D     = fracs.dist_mat(metric = metric)
    ## write distance matrix
    out_file = kwargs.get('out_file', 'sample_dist_' + metric +'.out')
    D.writetxt(out_file)    
    print 'wrote ' + out_file
    print 'Done!'
    

if __name__ == '__main__':
    ## parse input arguments
    from optparse import OptionParser
    kwargs = {}
    usage  = ('Compute the distance matrix between samples.\n' 
              'By default uses the the square-root of the Jensen-Shannon divergence.\n' 
              'distance matrix is written out as txt files. \n'
              'Requires the scipy.cluster.hierarchy module!\n'
              '\n'
              'Usage:   python SampleDist.py counts_file [options]\n'
              'Example: python SampleDist.py example/fake_data.txt -m JSsqrt -o my_dist_mat.out')
    parser = OptionParser(usage)
    parser.add_option("-m", "--metric", dest="metric", default='JSsqrt',
                      help="Distance metric to be utilized. JSsqrt (default) | any metric supported by scipy.cluster.hierarchy.")
    parser.add_option("-o", "--out_file", dest="out_file", type = 'str',
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs}, 
                      help="File to which distance matrix will be written.")
    (options, args) = parser.parse_args()
    counts_file     = args[0]
    metric          = options.metric
    ## write sample distance
    Run(counts_file, metric = metric, **kwargs)
    

