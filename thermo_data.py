
import sys, os

rmgpy_dir = os.environ.get('RMGpy')
if rmgpy_dir is None: raise RuntimeError("Please set RMGpy environment variable")
sys.path.append(os.path.join(rmgpy_dir,'source'))

import rmg
import rmg.thermo, rmg.structure, rmg.species

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
	os.path.exists('RMG_Database') or os.makedirs('RMG_Database')
	dictionary = file('Dictionary.txt','w')
	library = file('Library.txt','w')
	for s in list_of_species:
		dictentry, thermoentry = s.toRMGjava()
		dictionary.write(dictentry+"\n")
		library.write(thermoentry)
		print thermoentry
	dictionary.close()
	library.close()


# An example from cantera:
#species(name = "N2",
#    atoms = " N:2 ",
#    thermo = (
#       NASA( [  200.00,  1000.00], [  3.531005280E+00,  -1.236609880E-04, 
#               -5.029994330E-07,   2.435306120E-09,  -1.408812350E-12,
#               -1.046976280E+03,   2.967470380E+00] ),
#       NASA( [ 1000.00,  6000.00], [  2.952576370E+00,   1.396900400E-03, 
#               -4.926316030E-07,   7.860101950E-11,  -4.607552040E-15,
#               -9.239486880E+02,   5.871887620E+00] )
#             ),
#    note = "G 8/02"
#       )
#

franklins_method = """I did lots of cool stuff."""

acetylene = species(name = "acetylene",
    chemgraph = """
    1 C 0 {2,T} 
    2 C 0 {1,T} 
                """,
    thermo = (
       NASA( [  200.00,  1331.30], [  3.531005280E+00,  -1.236609880E-04, 
               -5.029994330E-07,   2.435306120E-09,  -1.408812350E-12,
               -1.046976280E+03,   2.967470380E+00] ),
       NASA( [ 1331.30,  6000.00], [  2.952576370E+00,   1.396900400E-03, 
               -4.926316030E-07,   7.860101950E-11,  -4.607552040E-15,
               -9.239486880E+02,   5.871887620E+00] )
             ),
    note = franklins_method + "Delta E down = 4"
       )


WriteRMGjava(_specieslist)

