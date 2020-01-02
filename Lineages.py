'''
Created on Jun 27, 2011

@author: jonathanfriedman
'''

import cPickle as pickle
import numpy as np

_levels         = ['k', 'p', 'c', 'o', 'f', 'g' , 's']
_levels_RDP     = {'domain':'k', 'phylum':'p', 'class':'c', 'order':'o', 'family':'f', 'genus':'g','species':'s'}
_unassigned_str = 'unassigned'

class Lineage(object):
    '''
    Class containing information regarding the taxonomic 
    assignment of an OTU.
    
    Class Attributes
    -----------------
    lin_str : str 
        The lineage line used to create self.
    lin : dict 
        Lineage assignments. keyed by taxonomic level 
        (1 letter abbreviations. See Lineages._levels)
    conf : dict 
        confidence level of the assignments. 
        Same keys as self.lin 
    '''

    def __init__(self, id=None, lin_str=None, format='QIIME'):
        '''        
        
        Parameters
        ----------
        id : hashable
            id of the current lineage. typically an OTU id.
        lin_str : str 
            line containing the lineage assignment information.
        format : {'QIIME' (default)| 'HMP' | 'RDP'} 
            The format of the lin_str.
            - QIIME : pairs of '[levelName]__[assignment]' 
              separated by semicolons. Level names are displayed even for 
              unassigned levels. Examples:
               * `k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__;s__`
               * `k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__;f__;g__;s__`         
            - HMP : pairs of '[assignment]([confidence])' separated by 
            semicolons. Unassigned levels are omitted. Examples:
               * `Root(100);Bacteria(100);"Firmicutes"(100);"Bacilli"(100);Bacillales(100);Bacillaceae(99);Bacillus(99);`
               * `Root(100);Bacteria(100);"Firmicutes"(100);"Clostridia"(100);Clostridiales(100);`  
            - RDP : triplets of '[assignment]\\t[levelName]\t[confidence]',
              separated by tabs. Unassigned levels are blanks, and all 
              tabs are retained! Examples:
               * `Bacteria\\tdomain\\t0.98\\tOD1\\tphylum\\t0.47\\t\\t\\t\\t\\t\\t\\t\\t\\t\\tOD1_genera_incertae_sedis\\tgenus\\t0.47`
        '''
        ## init id
        if not id: id = hash(self)
        self.id = id
        ##init taxonomy
        n_levels = len(_levels)
        lin_str = lin_str.replace('\t',';') 
        self.lin_str  = lin_str
        self.lin      = dict.fromkeys(_levels, _unassigned_str) # lineage information
        self.conf     = dict.fromkeys(_levels, _unassigned_str) # assignment confidence
        if lin_str: 
            if format.upper() == 'QIIME':
                fields = lin_str.strip().split(';')
                for field in fields:
                    field = field.strip()
                    if field.lower() == 'root':
                        continue
                    level, lin = field.split('__')
                    if lin == '': lin = _unassigned_str
                    self.lin[level] = lin.strip('"')
            elif format.upper() == 'HMP':
                fields = lin_str.strip().split(';')[1:-1] #ignore first assignment to root.
                n_fields = len(fields)
                fields += ['unassigned(0)'] * (n_levels-n_fields)
                for i,level in enumerate(_levels):
                    temp = fields[i].split('(')
                    assignment = temp[0].strip('"')
                    if assignment == 'unclassified':
                        assignment = _unassigned_str
                    self.lin[level]  = assignment
                    self.conf[level] = int(temp[-1][:-1]) 
            elif format.upper() == 'RDP':
                fields = lin_str.strip().split(';')
                for i,level in enumerate(_levels):
                    try:
                        lin  = fields[3*i]
                        conf = float(fields[3*i+2])
                        if not lin:
                            lin  = _unassigned_str
                            conf = 0                            
                    except:
                        lin  = _unassigned_str
                        conf = 0
                    self.lin[level]  = lin
                    self.conf[level] = conf                      
            else:
                raise ValueError("Unsuppoerted format: '%s'." %format)     
            
    def __repr__(self): return repr(self.lin_str)
    
    def __eq__(self, other):
        if type(self) != type(other):
            return False
        attrs = ('lin','lin_str','id','conf')
        return all( getattr(self,attr)==getattr(other,attr) for attr in attrs )
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def get_assignment(self, level=None, best=False):
        '''
        If best is True, return the most resolved lineage assignment 
        and its level.
        '''
        if level is None: # return most resolved assignment
            for l in _levels[::-1]:
                if self.lin[l] != _unassigned_str:
                    return l + '__' + self.lin[l]
            return l + '__' + self.lin[l]
        else: 
            if best: # return assignment at given level if assigned, o.w. return most resolved assignment.
                if self.lin[level] == _unassigned_str: 
                    return self.get_assignment(level=None,best=best)
                else: 
                    return self.get_assignment(level=level,best=False)
            else: # return assignment at given level (assigned or not)
                return level + '__' + self.lin[level]
        
        
class Lineages(dict):
    '''
    A dictionary of Lineage objects with dedicated methods.
    '''
    ## override some dict methods to support checking that all items are of type Lineage.    
    ## See: http://stackoverflow.com/questions/2060972/subclassing-python-dictionary-to-override-setitem
    def __init__(self, *args, **kwargs):
        '''
        A subcalss of :class:`dict` for holding :class:`Lineage` 
        objects.
        '''
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        if not isinstance(value, Lineage):
            t = type(value) 
            raise TypeError("Item corresponding to key '%s' is of type '%s instead of a Lineage object'." %(key, t)) 
        if key != value.id:
            raise ValueError("Item with ID '%s' has key '%s' instead of its own ID." %(value.id,key))    
        super(Lineages, self).__setitem__(key, value)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("update expected at most 1 arguments, got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]
    
    @classmethod
    def from_txt(cls, file, n_skip=0, format='QIIME', **kwargs):
        '''
        Create lineages object from txt file containing 2 columns: 
        id and lineage string.
        '''
        f     = open(file,'r')
        lines = f.readlines()
        delim_def = '\t'
        if format.upper() == 'RDP':
            delim_def = '\t\t\t\t\t'
        delim = kwargs.get('delim', delim_def)
        n = len(delim)
        d = {}
        for line in lines[n_skip:]:
            l = line.strip().replace('-','')
            i = l.index(delim)
            d[l[:i]] = l[i+n:]  
#        d     = dict([ line.strip().split(delim) for line in lines[n_skip:]  ])
        return cls.from_dict(d, format)
    
    
    def to_txt(self, file, sort_fun=str):
        '''
        Write all lineage strings to file, sorted by ids. 
        '''
        ids    = sorted(self.keys(), key=sort_fun)
        f      = open(file,'w')
        header = 'ID' +'\t' + 'Lineage' + '\n'
        f.write(header)
        for i in ids:
            line = str(i) +'\t' + self[i].lin_str + '\n'
            f.write(line)        
        f.close()
        
    
    def to_pickle(self,file):
        ''' pickles into file'''
        f=open(file,'w')
        pickle.dump(self,f)
        f.close()
    
    save = to_pickle
    
    @classmethod    
    def from_pickle(cls,file):
        ''' unpickles from file'''
        f=open(file,'r')
        temp=pickle.load(f)
        f.close()
        return temp
    
    @classmethod
    def from_dict(cls, d, format = 'QIIME'):
        '''
        Make a Lineages object from a dictionary whose keys are ids 
        and values are lineage strings.
        '''
        lins = Lineages()
        for id, s in d.iteritems():
            lin = Lineage(id = id, lin_str = s, format = format)
            lins[id] = lin
        return lins
          
    def get_assignments(self,level, best=False, ids='all'):
        '''
        Get the assignment of all lineages at given taxonomic level.
        
        Parameters
        ----------
        level: str 
            Desired taxonomic level.
        best: bool 
            If True, return best assignment if less resolved than required level. 
            Output in formal level_assignment
        ids : iterable/str  
            Ids for which to get taxonomic assignment.
        
        Returns
        ------- 
        A list of assignments keyed by ids of :class:`Lineage` objects.
        '''
        if level not in _levels: 
            raise ValueError("Level '%s' is not one of the allowed taxonomic leveles: 'k', 'p', 'c', 'o', 'f', 'g' , 's'." %level)
        if isinstance(ids, str):
            if ids=='all':
                a = [l.get_assignment(level, best=best) for l in self.itervalues()]
                return a
            else:
                ids = (ids,)
        a = [self[id].get_assignment(level, best=best) for id in ids]
        return a 
    
    
    def get_ids(self, level, assignment, best=True, complement = False):
        '''
        Get the ids of all items with a given taxonomy.
        If complement is True, get the ids of otus with assignment 
        different than given. 
        '''
        if level not in _levels: raise ValueError("Level '%s' is not one of the allowed taxonomic leveles: 'k', 'p', 'c', 'o', 'f', 'g' , 's'." %level)        
        if complement: ids = filter(lambda i: self[i].get_assignment(level=level, best=best) != assignment, self.iterkeys())
        else:          ids = filter(lambda i: self[i].get_assignment(level=level, best=best) == assignment, self.iterkeys())
        return ids
        
    def get_assigned(self, level):
        '''
        Get the ids of items that have an assigned taxonmy at given level.
        '''
        if level not in _levels: raise ValueError("Level '%s' is not one of the allowed taxonomic leveles: 'k', 'p', 'c', 'o', 'f', 'g' , 's'." %level)
        return self.get_ids(level, level + '.' + _unassigned_str, complement = True, best=False)


    def filter(self, ids):
        '''
        Return new instance with only given ids
        '''
        ids = filter(lambda i:i in self, ids)
        new_d = dict([ (id,self[id]) for id in ids ])
        new = Lineages(new_d)
        return new


if __name__ == '__main__':
    pass
#    test_Lineage()
    
    
            
