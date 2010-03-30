#!/usr/bin/env python
# encoding: utf-8
import sys, os

rmgpy_dir = os.environ.get('RMGpy')
if rmgpy_dir is None: raise RuntimeError("Please set RMGpy environment variable")
sys.path.append(os.path.join(rmgpy_dir,'source'))

import rmg
import rmg.thermo, rmg.structure, rmg.species, rmg.data, rmg.reaction
import rmg.log as logging


_rates = []


def loadKineticsDatabases(databasePath, only_families=False):
    """
    Create and load the kinetics databases (reaction families).
    If only_families is a list like ['H_Abstraction'] then only families in this
    list will be loaded.
    
    Currently incompatible with RMG(java) database syntax.
    """
    rmg.reaction.kineticsDatabase = rmg.reaction.ReactionFamilySet()
    rmg.reaction.kineticsDatabase.load(databasePath, only_families=only_families)
    

def removeCommentFromLine(line):
    """
    Remove a C++/Java style comment from a line of text.
    """
    index = line.find('//')
    if index >= 0:
        line = line[0:index]
    return line
    

def loadKineticsDatabases(databasePath, only_families=False):
    """
    Create and load the kinetics databases (reaction families).
    If only_families is a list like ['H_Abstraction'] then only families in this
    list will be loaded.
    """
    rmg.reaction.kineticsDatabase = rmg.reaction.ReactionFamilySet()
    rmg.reaction.kineticsDatabase.load(databasePath, only_families=only_families)
    db = rmg.reaction.kineticsDatabase
    return db
        


# Some Python code for parsing the kinetics library entries


class group:
    def __init__(self, chemgraph):
        s = rmg.structure.Structure()
        s.fromAdjacencyList(chemgraph.strip(),addH=False, withLabel=True)
        self._structure = s
        self.name = chemgraph.strip().split('\n')[0]
        self.node_name = None
        
    def __str__(self):
        return self._structure.toAdjacencyList(label=self.name)
    def __repr__(self):
        return 'group("""%s""")'%self.__str__()
    def toSource(self):
        return '"""%s"""'%self.__str__()
    source = property(toSource)
        
    def toRMGjava(self):
        if self.node_name is None:
            raise Exception("Location in family tree has not been found. Please call locateInFamily(reaction_family)")
        return self.node_name
    java = property(toRMGjava)
        
    def locateInFamily(self,family):
        """
        Locate this group in the family tree; set self.node_name to the node's name, and return it.
        Return None if not found.
        """
        s = self._structure
        node_name = family.descendTree(s,s.getLabeledAtoms())
        logging.debug("I think this %s belongs at %s"%(s.toAdjacencyList(label=group.name), node_name) )
        if node_name is not None: logging.debug( "and its ancestors are"+str(family.tree.ancestors(node_name)) )
        logging.debug("")
        
        if node_name != self.name:
            logging.warning("I think %s is actulaly %s"%(group.name, node_name))
        
        self.node_name = node_name
        return node_name

def indent(n,string):
    """Indent the string by n tabs (of 2 spaces each)"""
    out = ""
    for line in string.split('\n'):
        out+="    "*n + line + '\n'
    out = out[:-1] # remove the final '\n'
    return out

class rate:
    """
   
    ==================  ========================================
    Attribute           Description
    ==================  ========================================
    `group1`            The first group
    `group2`       
    `group3`      
    `reactant1`
    `reactant2` 
    `kf`                The forward rate constant, an :class:`Arrhenius` object
    `product1`          (maybe?)
    `product2` 
    `kr`                
    `rank`              A rough indicator of the quality of the data; use one of the constants (e.g. GUESS)
    `old_id`            The id number in the original database, if applicable
    `short_comment`     A one-line description of the data
    `long_comment`      A multi-line description of the data
    ==================  ========================================
   
    """

    def __init__(self, 
                 group1 = None,
                 group2 = None,
                 group3 = None,
                 reactant1 = None,
                 reactant2 = None,
                 kf = None,
                 rank = None,
                 old_id = None,
                 short_comment = None,
                 long_comment = None,
                 history = None
                ):
                
        if group1: self.group1 = group(group1)
        if group2: self.group2 = group(group2)
        if group3: self.group3 = group(group3)
        if reactant1: self.reactant1 = reactant(reactant1)
        if reactant2: self.reactant2 = reactant(reactant2)
        
        self.kf = kf
        self.rank = rank
        self.old_id = old_id
        self.short_comment = short_comment
        if short_comment.find('\n')>=0:
            logging.error("Line break found in short_comment: %s"%short_comment)
            raise Exception("No line breaks in short_comments please!")
        self.long_comment = long_comment
        
        global _rates
        _rates.append(self)
        
        
    def toRMGjava(self):
        out = ""
        for group in self.groups:
            out+="%-20s"%group.node_name
        out += self.kf.toRMGjava()
        out += '    '+self.short_comment
        return out
    java = property(toRMGjava)
        
    def toSource(self):
        out='rate(\n'
        # groups
        for i,group in enumerate(self.groups):
            out+='    group%d=""""\n%s           """",\n'%(i+1,group)
        out+='    kf = %r\n'%self.kf
        if self.rank:
            out+='    rank = %s,\n'%self.rank
        if self.old_id:
            out+='    old_id = "%s",\n'%self.old_id
        if self.short_comment:
            out+='    short_comment = "%s",\n'%self.short_comment
        if self.long_comment:
            out+='    long_comment = """%s\n                   """,\n'%self.long_comment.rstrip()
        out+=')'
        return out
    source = property(toSource)
        
    def getGroups(self):
        """Return a list of the groups: [group1, group2, (group3)]"""
        groups=[self.group1, self.group2]
        if hasattr(self,'group3'): groups.append(self.group3)
        return groups
    groups = property(getGroups)

class Parameter:
    """A parameter, with units, and uncertainties"""
    def __init__(self, *args):
        if args[0].__class__ is tuple: # unpack if double-packed
            args = args[0]
        nargs = len(args)
        self.value = args[0]
        self.units = None
        self.uncertainty_type = None
        self.uncertainty = None
        if nargs>1:
            self.units = args[1]
        if nargs>2:
            self.uncertainty_type = args[2]
            if nargs<4: raise Exception("Can't have uncertainty type without uncertainty or vice versa:"+str(args))
            self.uncertainty = args[3]
    def toSource(self):
        if self.uncertainty is not None:
            return '(%s, "%s", "%s", %s)'%(self.value, self.units, self.uncertainty_type, self.uncertainty)
        if self.units:
            return '(%s, "%s")'%(self.value, self.units)
        return '%s'%self.value
    source = property(toSource)
    def toRMGjava(self, unitfactor = 1.0):
        """Returns a tuple of strings (value, uncertainty), multiplied by unitfactor if set"""
        v = float(self.value) * unitfactor
        v = '%8.3g'%v
        if self.uncertainty is None:
            u = '%8s'%0
        elif self.uncertainty_type == '+-':
            u = float(self.uncertainty) * unitfactor
            u = '%8.3g'%u
        elif self.uncertainty_type == '*/':
            u = '*%.4g'%self.uncertainty
            u = '%8s'%u # make it 8 wide
        else: 
            raise Exception("Unknown uncertainty type %r in %r"%(self.uncertainty_type,self))
        return v,u
    java = property(toRMGjava)
        
    def __repr__(self):
        return "Parameter(%s)"%self.source.strip('()')
        
class Arrhenius:   
    """
    specify (A,n,Ea) or (A,n,alpha,E0)
    """
    def __init__(self, *args, **kw):
        self.A = None
        self.n = None
        self.Ea = None
        self.alpha = None
        self.E0 = None
        if kw.has_key('A'): self.A = Parameter(kw['A'])
        if kw.has_key('n'): self.n = Parameter(kw['n'])
        if kw.has_key('Ea'): self.Ea = Parameter(kw['Ea'])
        if kw.has_key('alpha'): self.alpha = Parameter(kw['alpha'])
        if kw.has_key('E0'): self.E0 = Parameter(kw['E0'])
        
    def __repr__(self):
        if self.Ea is not None:
            out = """Arrhenius( A=%s,
               n=%s,
               Ea=%s
              )"""%(self.A.source, self.n.source, self.Ea.source)
        if self.E0 is not None:
            out = """Arrhenius( A=%s,
               n=%s,
               alpha=%s
               E0=%s
              )"""%(self.A.source,self.n.source, self.alpha.source, self.E0.source)
        return out

    def toRMGjava(self):
        if self.alpha is None:
            alpha = Parameter(0)
            E0 = self.Ea
            logging.debug("Converting Arrhenius into Evans-Polyani")
        else:
            alpha = self.alpha
            E0 = self.E0
            
        if self.A.units == "cm^3/mol/s":
            Aunitscale = 1.0
        else:
            raise NotImplementedError("Unit conversion not yet implemented")
        Av,Au = self.A.toRMGjava(Aunitscale)
        # n and alpha have no units to convert 
        assert self.n.units is None
        nv,nu = self.n.toRMGjava(1)
        assert alpha.units is None
        alphav, alphau = alpha.toRMGjava(1)
        # Activation energy
        if E0.units == "kcal/mol":
            Eunitscale = 1.0
        else:
            raise NotImplementedError("Unit conversion not yet implemented")
        E0v, E0u = E0.toRMGjava(Eunitscale)
        
        out = ('%s '*8)%( Av, nv, alphav, E0v,
                          Au, nu, alphau, E0u )
        return out
    java = property(toRMGjava)

     
     
     
def WriteRMGjava(list_of_rates, family_name):
    folder = os.path.join('RMG_Database','kinetics_groups',family_name)
    os.path.exists(folder) or os.makedirs(folder)
    library = file(os.path.join(folder,'rateLibrary.txt'),'w')
    logging.verbose("Writing library to "+folder)
    for r in list_of_rates:
        line = r.toRMGjava()
        library.write(line+'\n')
        logging.verbose(line)
    library.close()
    

G3CANTHERM = 2.8 # for example
QUANTUM = 3
GUESS = 5

################################################################

# A sample kinetics library entry for H abstraction

family_name = 'H_Abstraction'
r=rate(
  group1 = """
    H2
    1 *1 H 0 {2,S}
    2 *2 H 0 {1,S}
    """,
  group2 = """
    C_rad/H/OneDeC
    1 *3 C             1 {2,S} {3,S} {4,S}
    2    H             0 {1,S}
    3    {Cd,Ct,Cb,CO} 0 {1,S}
    4    Cs            0 {1,S}
    """,
  kf = Arrhenius( A=(5.34E+11,"cm^3/mol/s", "*/", 1.3),
                  n=0, 
                  Ea=(24.0,"kcal/mol", "+-", 2.0)  # or alpha and E0
                 ),
  rank = QUANTUM,
  old_id = "140",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment.",
  long_comment = """
    Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
    """,
  history = [('2010-03-29',"Translated from old database",'rwest@mit.edu')]
  )
  
  
r = rate(
  group1 = 
"""
C/H3/Cs
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cs 0 {1,S}
""",
  group2 = 
"""
O_rad/NonDeO
1 *3 O 1 {2,S}
2    O 0 {1,S}
""",
  kf = Arrhenius(A=(2.52E+13,"cm^3/mol/s"),
                 n=0,
                 alpha=0,
                 E0=(20.435,"kcal/mol")),
  rank = 5,
  old_id = "148",
  short_comment = "Curran et al. [8] Rate expressions for H atom abstraction from fuels.",
  long_comment = 
"""
[8] Curran, H.J.; Gaffuri, P.; Pit z, W.J.; Westbrook, C.K. Combust. Flame 2002, 129, 253.
Rate expressions for H atom abstraction from fuels.

pg 257 A Comprehensive Modelling Study of iso-Octane Oxidation, Table 1. Radical:HO2, Site: primary (a)
Verified by Karma James
"""
)

r=rate(
  group1 = 
"""
C/H3/Cd
1 *1 C 0 {2,S}, {3,S}, {4,S}, {5,S}
2 *2 H 0 {1,S}
3 H 0 {1,S}
4 H 0 {1,S}
5 Cd 0 {1,S}
""",
  group2 = 
"""
C_rad/H2/Ct
1 *3 C 1 {2,S}, {3,S}, {4,S}
2 H 0 {1,S}
3 H 0 {1,S}
4 Ct 0 {1,S}
""",
  kf = Arrhenius(A=(7.58E+12,"cm^3/mol/s"),
                 n=0,
                 alpha=0,
                 E0=(18.2,"kcal/mol")),
  rank = 3,
  old_id = "125",
  short_comment = "Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.",
  long_comment = 
"""
Mark Saeys, CBS-QB3 calculations, without hindered rotor treatment. Rate expression per H atom.
"""
)


logging.initialize(20,'database.log')
# shut it up for loading the database!
logger = logging.getLogger()
old_level = logger.getEffectiveLevel()
logger.setLevel(40)
data_source = os.path.join('kinetics_groups_source')
db = loadKineticsDatabases( data_source )
logger.setLevel(old_level)


family_name = 'H abstraction'
family = db.families[family_name]

for r in _rates:
    for group in r.groups:
        group.locateInFamily(family)
        #group._structure.updateAtomTypes()

    logging.info("="*80)
    logging.info(r.source)
    logging.info(r.toRMGjava())

#import shutil
#shutil.copy(os.path.join(data_source,'kinetics_groups',family_name,'tree.txt'),os.path.join('RMG_Database','kinetics_groups',family_name,''))

WriteRMGjava(_rates, family_name)