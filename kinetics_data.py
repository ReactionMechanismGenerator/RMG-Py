#!/usr/bin/env python
# encoding: utf-8
import sys, os

rmgpy_dir = os.environ.get('RMGpy')
if rmgpy_dir is None: raise RuntimeError("Please set RMGpy environment variable")
sys.path.append(os.path.join(rmgpy_dir,'source'))


import rmg
import rmg.thermo, rmg.structure, rmg.species, rmg.data, rmg.reaction
import rmg.log as logging


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
		


_speciesnames = []
_specieslist = []

class structure():
	def __init__(self, chemgraph=None, SMILES=None):
		s = rmg.structure.Structure()
		if chemgraph is not None:
			s.fromAdjacencyList(chemgraph.strip(),addH=True, withLabel=False)
		if SMILES is not None:
			s.fromSMILES(SMILES)
		self._Structure = s
	def __repr__(self):
		return repr(self._Structure)
	def getFormula(self):
		return self._Structure.getFormula()
	def getSMILES(self):
		return self._Structure.toSMILES()
	def getChemgraph(self):
		return self._Structure.toAdjacencyList(strip_hydrogens=True)

## this creates the NASA class
from rmg.cantera_loader import NASA 
	# one of the methods is getRmgPolynomial()
		
class species():
	name=""
	structure=None
	note=""
	def __init__(self, 
				name,
				chemgraph=None,
				SMILES=None,
				thermo=None,
				note=None ):
		if name:
			global _speciesnames
			assert name not in _speciesnames, "Species %s already exists"%name
			self.name = name
		if chemgraph:
			self.structure = structure(chemgraph=chemgraph)
			for lin,lout in zip(chemgraph.split(),self.chemgraph.split()):
				if lin!=lout: 
					print "Warning - chemgraph has been changed from \n%s\n to\n%s\n"%(chemgraph,self.chemgraph)
					break
		if SMILES:
			self.structure = structure(SMILES=SMILES)
		if thermo:
			self.thermo = self.getRmgThermo(thermo)
		if note:
			self.note = note
			
		global _specieslist
		_specieslist.append(self)

	def getFormula(self):
		return self.structure.getFormula()
	formula = property(getFormula)
	def getSMILES(self):
		return self.structure.getSMILES()
	SMILES = property(getSMILES)
	def getChemgraph(self):
		return self.structure.getChemgraph()
	chemgraph = property(getChemgraph)
		
	def getRmgThermo(self,thermo):
		if all([isinstance(t, NASA) for t in thermo]): # all entries are NASA
			thermoData = rmg.thermo.ThermoNASAData(comment="Umm,")
			for thermoentry in thermo:
				poly = thermoentry.getRmgPolynomial()
				thermoData.addPolynomial(poly)
			return thermoData
		raise Exception("Can only deal with NASA thermo data")

	def toRMGjava(self):
		"""Get (dictionary_entry, library_entry) suitable for RMG Java"""
		dictentry = "%s\n"%self.name.strip()
		dictentry += self.chemgraph
		#//Name	H298	S298	Cp300	Cp400	Cp500	Cp600	Cp800	Cp1000	Cp1500	dH	dS	dCp	Comments
		J_to_cal = 0.239005736
		numbers = [self.thermo.getEnthalpy(298.15) * J_to_cal / 1000 ,
					self.thermo.getEntropy(298.15) * J_to_cal]
		numbers.extend([self.thermo.getHeatCapacity(T)*J_to_cal for T in [300, 400, 500, 600, 800, 1000, 1500]])
		numbers.extend([0, 0, 0])  # uncertainties!
		numbers_string = ("%8.3g"*12)%tuple(numbers)
		libraryentry = ("%-15s  %s %s")%(self.name, numbers_string, self.note)
		return dictentry, libraryentry
		
def WriteRMGjava(list_of_species):
	folder = os.path.join('RMG_Database','themo_libraries','latest')
	os.path.exists(folder) or os.makedirs(folder)
	dictionary = file(os.path.join(folder,'Dictionary.txt'),'w')
	library = file(os.path.join(folder,'Library.txt'),'w')
	for s in list_of_species:
		dictentry, thermoentry = s.toRMGjava()
		dictionary.write(dictentry+"\n")
		library.write(thermoentry)
		print thermoentry
	dictionary.close()
	library.close()

# Some Python code for parsing the kinetics library entries


class group:
    def __init__(self, chemgraph):
        s = rmg.structure.Structure()
        s.fromAdjacencyList(chemgraph.strip(),addH=False, withLabel=True)
        self._structure = s
        self.name = chemgraph.strip().split('\n')[0]
    def __str__(self):
        return self._structure.toAdjacencyList(label=self.name)
    def __repr__(self):
        return 'group("""%s""")'%self.__str__()
        
    def locateInFamily(self,family):
        """
        Locate this group in the family tree and set self.node_name to the node's name.
        """
        s = self._structure
        node_name = family.descendTree(s,s.getLabeledAtoms())
        logging.debug("I think this %s belongs at %s"%(s.toAdjacencyList(label=group.name), node_name) )
        if node_name is not None: logging.debug( "and its ancestors are"+str(fam.tree.ancestors(node_name)) )
        logging.debug("")
        
        if node_name != self.name:
            logging.warning("I think %s is actulaly %s"%(group.name, node_name))
        
        self.node_name = node_name
        return node_name
            
        
_rates = []
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
        if kf: self.kf = kf
        global _rates
        _rates.append(self)
    
    def getGroups(self):
        """Return a list of the groups: [group1, group2, (group3)]"""
        groups=[group1, group2]
        if group3: groups.append(group3)
        return groups
    groups = property(getGroups)

class Arrhenius:   
    """
    specify (A,n,Ea) or (A,n,alpha,E0)
    """
    def __init__(self, *args, **kw):
        if kw.has_key('A'): self.A = kw['A']
        if kw.has_key('n'): self.n = kw['n']





      
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
db = loadKineticsDatabases( os.path.join('kinetics_groups_source'))
logger.setLevel(old_level)


family = db.families['H abstraction']

for r in _rates:
    for group in [r.group1, r.group2]:
        group.locateInFamily(family)
        
    

