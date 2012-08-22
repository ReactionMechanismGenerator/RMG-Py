"""
authors: P Bhoorasingh, S Troiano
"""
import os
import logging
import numpy

from rmgpy.molecule import Molecule
from rmgpy.species import TransitionState
from molecule import QMMolecule, Geometry
from collections import defaultdict

import rdkit
from rdkit import DistanceGeometry
from rdkit.Chem.Pharm3D import EmbedLib

class QMReaction:
    
    file_store_path = 'QMfiles'
    if not os.path.exists(file_store_path):
        logging.info("Creating directory %s for mol files."%os.path.abspath(file_store_path))
        os.makedirs(file_store_path)
    
    def __init__(self):
        self.geometry = None
        self.transitionState = None
        
    def fixSortLabel(self, molecule):
        """
        This may not be required anymore. Was needed as when molecules were created, the
        rmg sorting labels would be set after where we tried to generate the TS.
        """
        sortLbl = 0
        for atom in molecule.atoms:
            atom.sortingLabel = sortLbl
            sortLbl += 1
        return molecule
    
    def atoms(self, mol):
        atoms = {}
        for atom in mol:
            args = atom.split()
            index = int(args.pop(0))
            if '*' in args[0]:
                label = args.pop(0)
            else:
                label = '  '
            type = args.pop(0)
            rad = args.pop(0)
            bonds = {}
            while args:
                bond = args.pop(0)[1:-1].split(',')
                bonds[int(bond[0])] = bond[1]
            atoms[index] = {'label': label, 'type': type,
                            'rad': rad, 'bonds': bonds}
        return atoms
    
    def adjlist(self, atoms):
        str = ''
        for key in atoms:
            atom = atoms[key]
            str += '\n{0:<{1}}{2}'.format(key,
                                          len('{0}'.format(max(atoms.keys()))) + 1,
                                          atom['label'])
            str += ' {0} {1}'.format(atom['type'], atom['rad'])
            for key0 in sorted(atom['bonds'].keys()):
                str += ' {' + '{0},{1}'.format(key0, atom['bonds'][key0]) + '}'
        return str.strip()
    
    def convert(self, rRMGMol, pRMGMol):
        rstrip = rRMGMol.toAdjacencyList().strip().splitlines()
        pstrip = pRMGMol.toAdjacencyList().strip().splitlines()
        reactant = self.atoms(rstrip)
        product = self.atoms(pstrip)
        output = {}
        swaps = {}
        import ipdb; ipdb.set_trace()
        for i in sorted(reactant.keys()):
            r_atom = reactant[i]
            if r_atom['label'].strip():
                for j in sorted(product.keys()):
                    p_atom = product[j]
                    if p_atom['label'] == r_atom['label']:
                        break
                output[i] = product.pop(j)
                if i != j:
                    swaps[j] = i
        for key in sorted(product.keys()):
            output[len(output) + 1] = product[key]
        for atom in output.values():
            bonds_in = {}
            bonds_out = {}
            for key in sorted(atom['bonds'].keys()):
                bonds_in[key] = atom['bonds'][key]
            for key in sorted(bonds_in.keys()):
                for old in swaps:
                    if key == old:
                        bonds_out[swaps[old]] = bonds_in[key]
                        bonds_in.pop(key)
                        break
            for key in sorted(bonds_in.keys()):
                bonds_out[key] = bonds_in[key]
            atom['bonds'] = {}
            for key in sorted(bonds_out.keys()):
                atom['bonds'][key] = bonds_out[key]
        return self.adjlist(output) + '\n'
    
    def getLabel(self, lbl1, lbl2, lbl3, lbl4):
        # Find the atom being transferred in the reaction
        if lbl1 == lbl3 or lbl1 == lbl4:
            atomLabel = lbl1
        elif lbl2 == lbl3 or lbl2 == lbl4:
            atomLabel = lbl2
            
        return atomLabel
    
    def getShiftedAtomLabel(self, actionList):
        for action in actionList:
            if action[0].lower() == 'form_bond':
                lbl1 = action[1]
                lbl2 = action[3]
            elif action[0].lower() == 'break_bond':
                lbl3 = action[1]
                lbl4 = action[3]
        
        return self.getLabel(lbl1, lbl2, lbl3, lbl4)
    
    def getSwitchedAtomLabel(self, actionList):
        for action in actionList:
            if action[0].lower() == 'form_bond':
                lbl1 = action[1]
                lbl2 = action[3]
            elif action[0].lower() == 'change_bond':
                lbl3 = action[1]
                lbl4 = action[3]
                
        return self.getLabel(lbl1, lbl2, lbl3, lbl4)
    
    def getBaseMolecules(self, atomLabel):
        try:
            self.reactants[0].getLabeledAtom(atomLabel)
            reactant1 = self.getDeepCopy(self.reactants[0])
            reactant2 = self.getDeepCopy(self.reactants[1])
        except ValueError:
            reactant1 = self.getDeepCopy(self.reactants[1])
            reactant2 = self.getDeepCopy(self.reactants[0])
            
        try:
            self.products[0].getLabeledAtom(atomLabel)
            product = self.getDeepCopy(self.products[0])
        except ValueError:
            product = self.getDeepCopy(self.products[1])
            
        return reactant1, reactant2, product
    
    def getGeometry(self, molecule):
        
        multiplicity = sum([i.radicalElectrons for i in molecule.atoms]) + 1
        geom = Geometry(molecule.toAugmentedInChIKey(), molecule, multiplicity)
        
        return geom, multiplicity
        
    def getRDKitMol(self, geometry):
        """
        Check there is no RDKit mol file already made. If so, use rdkit to make a rdmol from
        a mol file. If not, make rdmol from geometry.
        """ 
        if not os.path.exists(geometry.getCrudeMolFilePath()):
            geometry.generateRDKitGeometries()
        rdKitMol = rdkit.Chem.MolFromMolFile(geometry.getCrudeMolFilePath(), removeHs=False)      
        
        return rdKitMol
    
    def getDeepCopy(self, molecule):
        """
        Gets a deep copy of the molecule so that any edits to the molecule does not affect the
        original rmg molecule.
        """
        copyMol = molecule.copy(deep=True)
        
        return copyMol
    
    def chooseMol(self):
        if len(self.reactants) == 1:
            if self.reactants[0].atoms[0].sortingLabel == -1:
                self.reactants[0] = self.fixSortLabel(self.reactants[0])
            buildTS = self.getDeepCopy(self.reactants[0])
            actionList = self.family.forwardRecipe.actions
        else:
            if self.products[0].atoms[0].sortingLabel == -1:
                self.products[0] = self.fixSortLabel(reaction.products[0])
            buildTS = self.getDeepCopy(self.products[0])
            actionList = self.family.reverseRecipe.actions
        
        return buildTS, actionList
        
    def generateBoundsMatrix(self, molecule):
        """
        Uses rdkit to generate the bounds matrix of a rdkit molecule.
        """
        self.geometry, multiplicity = self.getGeometry(molecule)
        rdKitMol = self.getRDKitMol(self.geometry)
        boundsMatrix = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(rdKitMol)
        
        return rdKitMol, boundsMatrix, multiplicity
    
    def getBMParameters(self, reactant, product):
        for action in self.reaction.family.forwardRecipe.actions:
            if action[0].lower() == 'form_bond' or action[0].lower() == 'break_bond':
                atlbl1 = action[1]
                atlbl2 = action[3]
        
        if reactant.isCyclic():
            rdMol, bMatrix, mult = self.generateBoundsMatrix(reactant)
            lbl = (reactant.getLabeledAtom(atlbl1).sortingLabel, reactant.getLabeledAtom(atlbl2).sortingLabel)
            oRDMol, oBM, oMult = self.generateBoundsMatrix(product)
            other = (product.getLabeledAtom(atlbl1).sortingLabel, product.getLabeledAtom(atlbl2).sortingLabel)
            other = (oBM[other[0]][other[1]], oBM[other[1]][other[0]])
        elif product.isCyclic():
            rdMol, bMatrix, mult = self.generateBoundsMatrix(product)
            lbl = (product.getLabeledAtom(atlbl1).sortingLabel, product.getLabeledAtom(atlbl2).sortingLabel)
            oRDMol, oBM, oMult = self.generateBoundsMatrix(reactant)
            other = (reactant.getLabeledAtom(atlbl1).sortingLabel, reactant.getLabeledAtom(atlbl2).sortingLabel)
            other = (oBM[other[0]][other[1]], oBM[other[1]][other[0]])
        
        return rdMol, bMatrix, mult, lbl, other
        
    def editBondLimits(self, molecule, boundsMatrix, actionList):
        
        for action in actionList:
            lbl1 = action[1]
            atom1 = molecule.getLabeledAtom(lbl1)
            idx1 = atom1.sortingLabel
            if len(action) ==4:
                lbl2 = action[3]
                atom2 = molecule.getLabeledAtom(lbl2)
                idx2 = atom2.sortingLabel
            if action[0].lower() == 'change_bond':
                if action[2] == '1':
                    boundsMatrix[idx1][idx2] -= 0.15
                    boundsMatrix[idx2][idx1] -= 0.15
                elif action[2] == '-1':
                    boundsMatrix[idx1][idx2] += 0.15
                    boundsMatrix[idx2][idx1] += 0.15
            elif action[0].lower() == 'form_bond':
                boundsMatrix[idx1][idx2] -= 0.25
                boundsMatrix[idx2][idx1] -= 0.25
            elif action[0].lower() == 'break_bond':
                boundsMatrix[idx1][idx2] += 0.25
                boundsMatrix[idx2][idx1] += 0.25
            elif action[0].lower() == 'gain_radical':
                # may not be as simple as the others.
                pass
            elif action[0].lower() == 'lose_radical':
                pass
        
        # Smooth the bounds matrix to speed up the optimization
        rdkit.DistanceGeometry.DistGeom.DoTriangleSmoothing(boundsMatrix)
        
        return boundsMatrix
    
    def combineBoundsMatrices(self, boundsMat1, boundsMat2, atLbl1, atLbl2):
        # creates a bounds matrix for a TS by merging the matrices for its 2 reactants
        # Add bounds matrix 1 to corresponding place of the TS bounds matrix
        totSize = len(boundsMat1) + len(boundsMat2) - 1
        boundsMat = numpy.ones((totSize, totSize)) * 1000
        
        boundsMat[:len(boundsMat1),:len(boundsMat1)] = boundsMat1
        
        # Fill the bottom left of the bounds matrix with minima
        boundsMat[len(boundsMat1):, :len(boundsMat1)] = numpy.ones((len(boundsMat)-len(boundsMat1), len(boundsMat1))) * 1.07
        
        # Add bounds matrix 2, but it has to shift to the end of bounds matrix 1, and shift 
        # numbers for the reacting atom which has already been included from above
        boundsMat[len(boundsMat1):len(boundsMat1)+atLbl2, atLbl1] = boundsMat2[atLbl2, :atLbl2]
        boundsMat[atLbl1, len(boundsMat1):len(boundsMat1)+atLbl2] = boundsMat2[:atLbl2, atLbl2]
        boundsMat[atLbl1, len(boundsMat1)+atLbl2:] = boundsMat2[atLbl2, atLbl2+1:]
        boundsMat[len(boundsMat1)+atLbl2:, atLbl1] = boundsMat2[atLbl2+1:, atLbl2]
        
        # Remove all the parts of the transfered atom from the second bounds matrix
        # Incorporate the rest into the TS bounds matrix
        boundsMat2 = numpy.delete(numpy.delete(boundsMat2, atLbl2, 1), atLbl2, 0)
        boundsMat[-len(boundsMat2):, -len(boundsMat2):] = boundsMat2
        
        return boundsMat
    
    def matchAtoms(self, reactant, product):
        
        newadjlist = reactant.toAdjacencyList().strip().splitlines()
        radjlist = reactant.toAdjacencyList().strip().splitlines()
        padjlist = product.toAdjacencyList().strip().splitlines()
        rdict = {}
        for line in radjlist:
            if line.find('*') > -1:
                rdict[line.split()[1]] = int(line.split()[0])
        pdict = {}
        for line in padjlist:
            if line.find('*') > -1:
                pdict[line.split()[1]] = int(line.split()[0])
        
        rlabeldict = defaultdict(list) # label: (atom idx, atom symbol, bond type, radical)
        plabeldict = defaultdict(list)
        for label in sorted(reactant.getLabeledAtoms()):
            ridx = rdict[label] - 1
            pidx = pdict[label] - 1
            assert newadjlist[ridx].find('*') > -1
            line = newadjlist.pop(ridx)
            otherline = padjlist[pidx]
            for cut in otherline.split('{'):
                if cut.find('}') > -1:
                    check = int(cut.split(',')[0]) - 1
                    seg = padjlist[check]
                    if seg.find('*') > -1:
                        new = str(rdict[seg.split()[1]])
                        cut = cut.replace(cut.split(',')[0], new)
                        plabeldict[label].append('{' + cut)
            for cut in line.split('{'):
                if cut.find('}') > -1:
                    check = int(cut.split(',')[0]) - 1
                    seg = radjlist[check]
                    if seg.find('*') > -1:
                        rlabeldict[label].append('{' + cut)
            for bond in rlabeldict[label]:
                line = line.replace(bond, '')
            for bond in plabeldict[label]:
                line = line + ' ' + bond.replace(' ', '')
            newadjlist.insert(ridx, line)
        
        return newadjlist
        
    def twoEnded(self):
        fam = self.reaction.family.name
        if len(self.reactants) == len(self.products):
            if len(self.reactants) == 1:
                reactant = self.getDeepCopy(self.reactants[0])
                product = self.getDeepCopy(self.products[0])
            else:
                reactant1 = self.getDeepCopy(self.reactants[0])
                reactant2 = self.getDeepCopy(self.reactants[1])
                product1 = self.getDeepCopy(self.products[0])
                product2 = self.getDeepCopy(self.products[1])
                reactant = reactant1.merge(reactant2)
                product = product1.merge(product2)
        else:
            if len(self.reactants) == 1:
                reactant = self.getDeepCopy(self.reactants[0])
                product1 = self.getDeepCopy(self.products[0])
                product2 = self.getDeepCopy(self.products[1])
                product = product1.merge(product2)
            else:
                reactant1 = self.getDeepCopy(self.reactants[0])
                reactant2 = self.getDeepCopy(self.reactants[1])
                product = self.getDeepCopy(self.products[0])
                reactant = reactant1.merge(reactant2)
        
        newadjlist = self.matchAtoms(reactant, product)
        newadjlist = self.adjlist(self.atoms(newadjlist))
        test = Molecule()
        if product.isIsomorphic(test.fromAdjacencyList(newadjlist)):
            product = product.fromAdjacencyList(newadjlist)
        else:
            logging.info("Couldn't generate transition state geometry")
            
        return reactant, product
    
    def stretchBond(self, boundsMatrix, lbl):
        boundsMatrix[lbl[0]][lbl[1]] += 0.24
        boundsMatrix[lbl[1]][lbl[0]] += 0.24
        
        DistanceGeometry.DistGeom.DoTriangleSmoothing(boundsMatrix)
        
        return boundsMatrix
    
    def getTSBMatrix(self):
        
        reactant = self.getDeepCopy(self.reactants[0])
        product = self.getDeepCopy(self.products[0])
        if reactant.atoms[0].sortingLabel == -1:
            reactant = self.fixSortLabel(reactant)
        if product.atoms[0].sortingLabel == -1:
            product = self.fixSortLabel(product)
            
        rdMol, tsBM, mult, lbl, other = self.getBMParameters(reactant, product)
        
        return rdMol, tsBM, mult, lbl, other
        
    def generateGeometry(self):
        """
        This is for reactions in the family of intra R Add Exo/Endo-cyclic
        """
        # A --> B or A + B --> C + D
        if len(self.reactants) == len(self.products):
            actionList = self.family.forwardRecipe.actions
            # 1 reactant
            if len(self.reactants) == 1:
                reactant = self.getDeepCopy(self.reactants[0])
                product = self.getDeepCopy(self.products[0])
                adjlist = self.convert(reactant, product)
                if product == product.fromAdjacencyList(adjlist):
                    product = product.fromAdjacencyList(adjlist)
                else:
                    logging.info("Couldn't generate transition state geometry")
                    
                reactantRDMol, boundsMatR, multiplicityR = self.generateBoundsMatrix(reactant)
                productRDMmol, boundsMatP, multiplicityP = self.generateBoundsMatrix(product)
                
                # Edit one bounds matrix with values from the other
                atLbl = self.getSwitchedAtomLabel(actionList)
                
                rAtLbl = reactant.getLabeledAtom(atLbl).sortingLabel
                pAtLbl = product.getLabeledAtom(atLbl).sortingLabel
                # to do this properly I will have to match atoms from the product and reactant
                boundsMatR[:][rAtLbl] = boundsMatP[:][pAtLbl]
                boundsMatR[rAtLbl][:] = boundsMatP[pAtLbl][:]
                
            # 2 reactants
            else:
                atLbl = self.getShiftedAtomLabel(actionList)
                
                # Derive the bounds matrix from the reactants and products
                reactant, reactant2, product = self.getBaseMolecules(atLbl)
                
                # Check for sorting labels
                if reactant.atoms[0].sortingLabel != 0:
                    reactant = self.fixSortLabel(reactant)
                if product.atoms[0].sortingLabel != 0:
                    product = self.fixSortLabel(product)
                
                # Generate the bounds matrices for the reactant and product with the transfered atom
                reactantRDMol, boundsMatR, multiplicityR = self.generateBoundsMatrix(reactant)
                productRDMmol, boundsMatP, multiplicityP = self.generateBoundsMatrix(product)
                
                rAtLbl = reactant.getLabeledAtom(atLbl).sortingLabel
                pAtLbl = product.getLabeledAtom(atLbl).sortingLabel
                
                # Get the total size of the TS bounds matrix and initialize it
                boundsMat = self.combineBoundsMatrices(boundsMatR, boundsMatP, rAtLbl, pAtLbl)
                
                # Merge the reactants to generate the TS template
                buildTS = reactant.merge(reactant2)
                self.fixSortLabel(buildTS)
                boundsMat = self.editBondLimits(buildTS, boundsMat, actionList)
                
                self.geometry, multiplicity = self.getGeometry(buildTS)
                
                self.rdmol = self.getRDKitMol(self.geometry)
                
        # A --> B + C or A + B --> C
        else:
            # The single species is used as the base for the transition state 
            buildTS, actionList = self.chooseMol()
            
            # Generate the RDKit::Mol from the RMG molecule and get the bounds matrix
            self.rdmol, boundsMat, multiplicity = self.generateBoundsMatrix(buildTS)
            
            # Alter the bounds matrix based on the reaction recipe
            boundsMat = self.editBondLimits(buildTS, boundsMat, actionList)
             
        # Optimize the TS geometry in place, outputing the initial and final energies
        try:
            rdkit.Chem.Pharm3D.EmbedLib.OptimizeMol(self.rdmol, boundsMat, maxPasses = 10)
        except RuntimeError:
            pass
        
        return buildTS, boundsMat, multiplicity