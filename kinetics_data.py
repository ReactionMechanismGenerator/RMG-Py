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
    
    db = rmg.reaction.ReactionFamilySet()
    rmg.reaction.kineticsDatabase = db
    
    datapath = os.path.abspath(databasePath)
    logging.info('Loading reaction family databases from %s.' % datapath)
    
    # Load the families from kinetics/families.txt
    familyList = []
    ffam = open(os.path.join(datapath,'kinetics_groups','families.txt'), 'r')
    for line in ffam:
        line = rmg.data.removeCommentFromLine(line).strip()
        if len(line) > 0:
            items = line.split()
            items[0] = int(items[0])
            familyList.append(items)
    ffam.close()
    
    # loop through the families
    families = {}
    for index, status, label in familyList:
        path = os.path.join(datapath, 'kinetics_groups', label)
        if os.path.isdir(path): # and status.lower() == 'on':
            # skip families not in only_families, if it's set
            if only_families and label not in only_families: continue
            
            logging.verbose('Loading reaction family %s from %s...' % (label, datapath))
            family = rmg.reaction.ReactionFamily(label)
            #family.load(path)  # this would load the whole family. We just want the dictionary and tree
            
            # Generate paths to files in the database
            dictstr = os.path.join(path, 'dictionary.txt')
            treestr = os.path.join(path, 'tree.txt')
            #: The path of the database that was loaded.
            family._path = path
            
            rmg.data.Database.load(family, dictstr, treestr, '')
            
            # Check for well-formedness
            if not family.isWellFormed():
                raise data.InvalidDatabaseException('Database at "%s" is not well-formed.' % (path))                
                    
            db.families[family.label] = family

    # rmg.reaction.kineticsDatabase.load(databasePath, only_families=only_families)
    db = rmg.reaction.kineticsDatabase
    return db
        


# Some Python code for parsing the kinetics library entries


class Group:
    def __init__(self, chemgraph):
        chemgraph = chemgraph.strip()
        self.name = chemgraph.splitlines()[0]
        self.node_name = None
        if not chemgraph.splitlines()[1].strip().startswith('1'):
            #logging.debug("I think this is a compound structure: %s"%chemgraph)
            self.compound_structure=True
            self.chemgraph = chemgraph
            return
        else:
            self.compound_structure=False
        s = rmg.structure.Structure()
        s.fromAdjacencyList(chemgraph.strip(),addH=False, withLabel=True)
        self._structure = s
        
    def __str__(self):
        if self.compound_structure:
            return self.chemgraph+'\n'
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
        
        if self.compound_structure:
            if family.dictionary.has_key(self.name):
                logging.info("%s is a compound group so I'm just using its name and trusting that's right."%self.name)
                self.node_name = self.name
            else:
                logging.warning("%s is a compound group but there's no node of the same name, so I'm setting it to None!"%self.name)
                self.node_name = None
            return self.node_name
        if self.name in ['R6_SSSDS']:
            logging.debug("I'm about to try and match "+self.name)
        s = self._structure
        node_name = family.descendTree(s,s.getLabeledAtoms())
        logging.debug("I think this %s belongs at %s"%(s.toAdjacencyList(label=self.name), node_name) )
        if node_name is not None: logging.debug( "and its ancestors are"+str(family.tree.ancestors(node_name)) )
        logging.debug("")
        
        if node_name != self.name:
            logging.warning("WARNING! in %s family I think %s is actually at node %s"%(family.label, self.name, node_name))
            if node_name:
                logging.verbose(node_name + ":")
                if family.dictionary[node_name].__class__ is rmg.structure.Structure:
                    logging.verbose(family.dictionary[node_name].toAdjacencyList())
                else: 
                    logging.verbose(str(family.dictionary[node_name]))
                #import pdb; pdb.set_trace()
            if family.dictionary.has_key(self.name):
                logging.warning("There's a node called '%s' so I'm using that."%self.name)
                node_name = self.name
                if family.tree.parent.has_key(node_name):
                    logging.warning( "Its ancestors are"+str(family.tree.ancestors(node_name)) )
                else:
                    logging.warning("It is not in the tree!!!")
                logging.verbose(self.name+": \n"+family.dictionary[self.name].toAdjacencyList() )

        self.node_name = node_name
        return node_name
        
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
        if self.uncertainty:
            return '(%g, "%s", "%s", %g)'%(self.value, self.units, self.uncertainty_type, self.uncertainty)
        if self.units:
            return '(%g, "%s")'%(self.value, self.units)
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
        A_val,A_unc = self.A.toRMGjava(Aunitscale)
        # n and alpha have no units to convert 
        assert self.n.units is None
        n_val,n_unc = self.n.toRMGjava(1)
        assert alpha.units is None
        alpha_val, alpha_unc = alpha.toRMGjava(1)
        # Activation energy
        if E0.units == "kcal/mol":
            Eunitscale = 1.0
        else:
            raise NotImplementedError("Unit conversion not yet implemented")
        E0_val, E0_unc = E0.toRMGjava(Eunitscale)
        
        out = ("%8s "*8)%(
                    A_val, n_val, alpha_val, E0_val, A_unc, n_unc, alpha_unc, E0_unc )
        return out
    java = property(toRMGjava)
    header = ("%8s "*8)%('A','n','alpha','E0','DA','Dn','Dalpha','DE0')
    

def indent(n,string):
    """Indent the string by n tabs (of 2 spaces each)"""
    out = ""
    for line in string.split('\n'):
        out+="    "*n + line + '\n'
    out = out[:-1] # remove the final '\n'
    return out

class rate:
    """
    Class used to define reaction rates in the kinetics library.
   
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
                 temperature_range = None,
                 rank = None,
                 old_id = None,
                 short_comment = None,
                 long_comment = None,
                 history = None
                ):
                
        if group1: self.group1 = Group(group1)
        if group2: self.group2 = Group(group2)
        if group3: self.group3 = Group(group3)
        if reactant1: self.reactant1 = reactant(reactant1)
        if reactant2: self.reactant2 = reactant(reactant2)
        
        self.kf = kf
        self.temperature_range = temperature_range
        self.rank = rank
        self.old_id = old_id
        self.short_comment = short_comment
        if short_comment.find('\n')>=0:
            logging.error("Line break found in short_comment: %s"%short_comment)
            raise Exception("No line breaks in short_comments please!")
        self.history = history
        self.long_comment = long_comment
        
        global _rates
        _rates.append(self)
        
        
    def toRMGjava(self):
        out = ""
        for group in self.groups:
            out+="%-20s "%group.node_name
        if self.temperature_range[1] is None:
            out += "  %6d    "%self.temperature_range[0]
        else:
            out += " %4d-%-5d "%self.temperature_range
            
        out += self.kf.toRMGjava()
        out += " %d "%int(self.rank or 0)
        out += '    '+self.short_comment
        return out
    java = property(toRMGjava)
    header =( "//#  %-18s "%'Group1' + 
              "%-20s "%'Group2' +
              " TempRange  " +
              Arrhenius.header +
              " Rank " +
              " Comment" )
        
    def toSource(self):
        out='rate(\n'
        # groups
        for i,group in enumerate(self.groups):
            out+='    group%d =\n"""\n%s""",\n'%(i+1,group)
        out+='    kf = %r\n'%self.kf
        if self.temperature_range:
            out+='    temperature_range = %s,\n'%repr(self.temperature_range)
        if self.rank is not None:
            out+='    rank = %r,\n'%self.rank
        if self.old_id:
            out+='    old_id = %r,\n'%self.old_id
        if self.short_comment:
            out+='    short_comment = %r,\n'%self.short_comment
        if self.long_comment:
            out+='    long_comment =\n"""\n%s\n""",\n'%self.long_comment.strip()
        if self.history:
            out+='    history = %s\n'%str(self.history)
        out+=')'
        return out
    source = property(toSource)
        
    def getGroups(self):
        """Return a list of the groups: [group1, group2, group3] or as few as exist"""
        groups=[self.group1]
        if hasattr(self,'group2'): groups.append(self.group2)
        if hasattr(self,'group3'): groups.append(self.group3)
        return groups
    groups = property(getGroups)




     
def WriteRMGjava(list_of_rates, output_folder, unread_lines=None, header=None ):
    os.path.exists(output_folder) or os.makedirs(output_folder)
    library = file(os.path.join(output_folder,'rateLibrary.txt'),'w')
    logging.verbose("Writing library to "+output_folder)
    if unread_lines:
        library.write(unread_lines+'\n')
    if header:
        library.write(header+'\n')
    for index,r in enumerate(list_of_rates):
        line = '%-4d '%index
        line += r.toRMGjava()
        library.write(line+'\n')
        logging.verbose(line)
    library.close()
    
def WriteSource(list_of_rates, output_folder, family_name=None, file_header='', unread_lines=None):
    os.path.exists(output_folder) or os.makedirs(output_folder)
    library = file(os.path.join(output_folder,'library.py'),'w')
    logging.verbose("Writing library to "+output_folder)
    
    library.write('# encoding: utf-8\n')
    library.write('header = """\n%s\n"""\n\n'%file_header.replace('"','\\"').strip())
    
    library.write('reaction_family_name = "%s"\n'%family_name)
    
    library.write('# These lines were in the RMG library file but were not translated into anything useful:\n')
    library.write('unread_lines = ')
    if unread_lines:
        library.write('"""'+unread_lines.replace('"','\\"')+'"""\n')
    else:
        library.write('None\n')
    
    
    for index,r in enumerate(list_of_rates):
        library.write("# Number %d\n"%index)
        library.write(r.source+'\n\n')
        #logging.verbose(line)
        
    library.write('\n')
    library.close()
    

G3CANTHERM = 2.8 # for example
QUANTUM = 3
GUESS = 5

################################################################

def main():
    pass
    

import shutil
logging.initialize(5,'kinetics_data.log')  # <<< level in general

logger = logging.getLogger()
old_level = logger.getEffectiveLevel()
logger.setLevel(15)  # <<< level while reading database 
data_source = os.path.join('input','RMG_database')
db = loadKineticsDatabases( data_source )
logger.setLevel(old_level)

for family_name,family in db.families.iteritems():
    logging.info("Now processing family: %s"%family_name)
    _rates = []
    
    library_file_name = os.path.join(data_source,'kinetics_groups',family_name,'library.py')
    library_file = open(library_file_name)
    context = { '_rates': _rates,
              'rate': rate,
              'Arrhenius':Arrhenius,
              'Parameter':Parameter,
              "__builtins__":None
              }
    logging.info("Importing %s"%library_file_name)
    # TODO: check library_file is safe. Perhaps with http://code.activestate.com/recipes/496746-restricted-safe-eval/
    # although I think limiting the context may be ample protection:
    try:
        exec library_file in {"__builtins__":None}, context
    except NameError, e:
        logging.error("Looks like the file %s had an illegal operator in it:"%library_file_name)
        logging.error(e)
        raise
    finally:
        library_file.close()
    _rates = context['_rates']
    unread_lines = context['unread_lines']
    file_header = context['header']
    
    assert family_name==context['reaction_family_name']
    
    for r in _rates:
        #logging.verbose("="*80)
        for group in r.groups:
            group.locateInFamily(family)
            #group._structure.updateAtomTypes()
            
        #logging.verbose(r.source)
        #logging.verbose('#RMG-Java: '+r.toRMGjava())
        
    output_folder = os.path.join('output','RMG_Database','kinetics_groups',family_name)
    os.path.exists(output_folder) or os.makedirs(output_folder)
    WriteRMGjava(_rates, output_folder, unread_lines=unread_lines, header = rate.header)
    WriteSource(_rates, output_folder, file_header=file_header, unread_lines=unread_lines, family_name=family_name)
    
    for file_to_copy in ['tree.txt', 'dictionary.txt', 'reactionAdjList.txt']:
        logging.debug("Copying %s into destination folder"%file_to_copy)
        shutil.copy( os.path.join(data_source,'kinetics_groups',family_name,file_to_copy),
                    output_folder)
