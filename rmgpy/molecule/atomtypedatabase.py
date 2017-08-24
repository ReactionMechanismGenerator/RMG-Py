'''
Created on Sep 4, 2014

@author: nickvandewiele
'''
class AbstractAtomType(object):
    def __init__(self, element = None, label=None, double=-1, triple=-1, benzene=-1, lp=-1):
        self.element = element
        self.label = label
        self.double = double
        self.triple = triple
        self.benzene = benzene
        self.lp = lp

class Column4(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.lp = 0
            
class Column5(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.lp = 1

class Column6(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.lp = 2        

class Xs(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.double, self.triple, self.benzene = 0, 0, 0
            self.label = 's'
            
class Xd(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.double, self.triple, self.benzene = 1, 0, 0
            self.label = 'd'
            
class Xdd(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.double, self.triple, self.benzene = 2, 0, 0
            self.label = 'dd'
            
class Xt(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.double, self.triple, self.benzene = 0, 1, 0
            self.label = 't'
                        
class Xb(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.double, self.triple, self.benzene = 0, 0, 2
            self.label = 'b'            

class Xbf(AbstractAtomType):
    def __init__(self, *args, **kwargs):
            super(self.__class__, self).__init__(*args, **kwargs)
            self.double, self.triple, self.benzene = 0, 0, 3
            self.label = 'bf'
                            
def create_atom_types():
    atomtypes = []
    
    #tetravalent:
    
    tetravalent = []
    for type in [Xs, Xd, Xdd, Xt, Xb, Xbf]:
        tetravalent.extend(create_types(type, ['C', 'Si']))
    
    for at in tetravalent: at.lp = 0
    
    atomtypes.extend(tetravalent)
            
    #bivalent:
    bivalent = []
    for type in [Xs, Xd]:
        #bivalent.extend(create_types(type, ['O', 'S']))
        bivalent.extend(create_types(type, ['O']))
    
    for at in bivalent: at.lp = 2
    
    atomtypes.extend(bivalent)
    
    #trivalent nitrogen:
    trivalent_N = []
    for type in [Xs, Xd, Xt, Xb]:
        trivalent_N.extend(create_types(type, ['N'], ['N3']))
    
    for at in trivalent_N: at.lp = 1
    atomtypes.extend(trivalent_N)
    
    #pentavalent nitrogen:
    pentavalent_N = []
    for type in [Xs, Xd, Xdd, Xt, Xb]:
        pentavalent_N.extend(create_types(type, ['N'], ['N5']))
    
    for at in pentavalent_N: at.lp = 0
    atomtypes.extend(pentavalent_N)
    
    return atomtypes
    
def create_types(Type, elements, labels=None):
    if labels is None:
        labels = elements
    atomtypes = []
    for el, label in zip(elements, labels):
        at = Type(element=el)
        at.label = label+at.label
        atomtypes.append(at)
    return atomtypes
