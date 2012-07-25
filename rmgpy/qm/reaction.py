def generateTransitionState(self, reaction):
    """
    make TS geometry
    """
    from rmgpy.molecule import Molecule
    from rmgpy.species import TransitionState
    import rdkit
    from rdkit.Chem.Pharm3D import EmbedLib
    
    def fixSortLabel(molecule):
        sortLbl = 0
        for atom in molecule.atoms:
            atom.sortingLabel = sortLbl
            sortLbl += 1
        return molecule
    
    # A --> B or A + B --> C + D               
    if len(reaction.reactants) == len(reaction.products):
        # 1 reactant
        if len(reaction.reactants) == 1:
            pass
        
        # 2 reactants
        else:
            actionList = reaction.family.forwardRecipe.actions
            for action in actionList:
                if action[0].lower() == 'form_bond':
                    lbl1 = action[1]
                    lbl2 = action[3]
                elif action[0].lower() == 'break_bond':
                    lbl3 = action[1]
                    lbl4 = action[3]
            
            # Find the atom being transferred in the reaction
            if lbl1 == lbl3 or lbl1 == lbl4:
                lblAt = lbl1
            else:
                lblAt = lbl2
            
            # Derive the bounds matrix from the reactants and products
            try:
                reaction.reactants[0].getLabeledAtom(lblAt)
                reactant = reaction.reactants[0]
                reactant2 = reaction.reactants[1]
            except ValueError:
                reactant = reaction.reactants[1]
                reactant2 = reaction.reactants[0]
                
            try:
                reaction.products[0].getLabeledAtom(lblAt)
                product = reaction.products[0]
            except ValueError:
                product = reaction.products[1]
            
            # Merge the reactants to generate the TS template
            import ipdb; ipdb.set_trace()
            buildTS = reactant.merge(reactant2)
            
            # Check for rdkit molecules
            if reactant.rdMol == None:
                Molecule.generate3dGeometry(reactant)
            if product.rdMol == None:
                Molecule.generate3dGeometry(product)
            if buildTS.rdMol == None:
                Molecule.generate3dGeometry(buildTS)
                
            # Check for sorting labels
            if reactant.atoms[0].sortingLabel != 0:
                reactant = fixSortLabel(reactant)
            if product.atoms[0].sortingLabel != 0:
                product = fixSortLabel(product)
            
            # Generate the bounds matrices for the reactant and product with the transfered atom
            boundsMat1 = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(reactant.rdMol)
            boundsMat2 = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(product.rdMol)
            
            # Get the total size of the TS bounds matrix and initialize it
            totSize = len(boundsMat1) + len(boundsMat2) - 1
            boundsMat = numpy.ones((totSize, totSize)) * 1000
            
            # Add bounds matrix 1 to corresponding place of the TS bounds matrix
            boundsMat[:len(boundsMat1),:len(boundsMat1)] = boundsMat1
            
            # Fill the bottom left of the bounds matrix with minima
            boundsMat[len(boundsMat1):, :len(boundsMat1)] = numpy.ones((len(boundsMat)-len(boundsMat1), len(boundsMat1))) * 1.07
            
            # Add bounds matrix 2, but it has to shift to the end of bounds matrix 1, and shift 
            # numbers for the reacting atom which has already been included from above
            rAtLbl = reactant.getLabeledAtom(lblAt).sortingLabel
            pAtLbl = product.getLabeledAtom(lblAt).sortingLabel
            boundsMat[len(boundsMat1):len(boundsMat1)+pAtLbl, rAtLbl] = boundsMat2[pAtLbl, :pAtLbl]
            boundsMat[rAtLbl, len(boundsMat1):len(boundsMat1)+pAtLbl] = boundsMat2[:pAtLbl, pAtLbl]
            boundsMat[rAtLbl, len(boundsMat1)+pAtLbl+1:] = boundsMat2[pAtLbl, pAtLbl+1:]
            boundsMat[len(boundsMat1)+pAtLbl+1:, rAtLbl] = boundsMat2[pAtLbl+1:, pAtLbl]
            
            # Remove all the parts of the transfered atom from the second bounds matrix
            # Incorporate the rest into the TS bounds matrix
            boundsMat2 = numpy.delete(numpy.delete(boundsMat2, pAtLbl, 1), pAtLbl, 0)
            boundsMat[-len(boundsMat2):, -len(boundsMat2):] = boundsMat2
            
            #********what now!!!??!?
            
    # A --> B + C or A + B --> C
    else:
        # Fix the sorting label for the molecule if it has not been done.
        # Set the action list to forward or reverse depending on the species
        # the transition state is being built from.
        if len(reaction.reactants) == 1:
            if reaction.reactants[0].atoms[0].sortingLabel == -1:
                reaction.reactants[0] = fixSortLabel(reaction.reactants[0])
            buildTS = reaction.reactants[0]
            actionList = reaction.family.forwardRecipe.actions
        else:
            if reaction.products[0].atoms[0].sortingLabel == -1:
                reaction.products[0] = fixSortLabel(reaction.products[0])
            buildTS = reaction.products[0]
            actionList = reaction.family.reverseRecipe.actions
        
        # Generate the RDKit::Mol from the RMG molecule and get the bounds matrix
        if buildTS.rdMol == None:
            Molecule.generate3dGeometry(buildTS)
        boundsMat = rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(buildTS.rdMol)
        
        # Alter the bounds matrix based on the reaction recipe
        for action in actionList:
            lbl1 = action[1]
            atom1 = buildTS.getLabeledAtom(lbl1)
            idx1 = atom1.sortingLabel
            if len(action) ==4:
                lbl2 = action[3]
                atom2 = buildTS.getLabeledAtom(lbl2)
                idx2 = atom2.sortingLabel
            if action[0].lower() == 'change_bond':
                # bond was added
                if action[2] == '1':
                    # make the bond shorter
                    boundsMat[idx1][idx2] -= 0.15
                    boundsMat[idx2][idx1] -= 0.15
                # bond was removed
                elif action[2] == '-1':
                    # make bond shorter
                    boundsMat[idx1][idx2] += 0.15
                    boundsMat[idx2][idx1] += 0.15
            elif action[0].lower() == 'form_bond':
                # move them further
                boundsMat[idx1][idx2] -= 0.25
                boundsMat[idx2][idx1] -= 0.25
            elif action[0].lower() == 'break_bond':
                # move them closer
                boundsMat[idx1][idx2] += 0.25
                boundsMat[idx2][idx1] += 0.25
            elif action[0].lower() == 'gain_radical':
                pass
            elif action[0].lower() == 'lose_radical':
                pass
        
        # Keep the most stable conformer, remove the rest
        for conf in range(0, buildTS.rdMol.GetNumConformers()):
            if conf != buildTS.rdMolConfId:
                buildTS.rdMol.RemoveConformer(conf)
        
        # Smooth the bounds matrix to speed up the optimization
        # Optimize the TS geometry in place, outputing the initial and final energies
        rdkit.DistanceGeometry.DistGeom.DoTriangleSmoothing(boundsMat)
        try:
            rdkit.Chem.Pharm3D.EmbedLib.OptimizeMol(buildTS.rdMol, boundsMat, maxPasses = 10)
        except RuntimeError:
            pass
        
        # Ensure the spin multiplicity is ok for gaussian
        spinMult = 1
        for atom in buildTS.atoms:
            if atom.spinMultiplicity !=1:
                spinMult += 1
        
        # Need to output gaussian script
        title = reaction.reactants[0].toInChI().replace('/','_').strip('InChI=1S_') + reaction.reactants[1].toInChI().replace('/','_').strip('InChI=1S_')
        filename = title + '.gif'
        newPath = os.path.join(os.curdir, 'tsScripts')
        # Need to makedirs if gaussScripts not there
        fout = open(os.path.join(newPath, filename), "w")
        fout.write("# PM3 Opt=(TS, CalcAll, NoEigenTest)\n\n")
        fout.write(title+"\n\n")
        fout.write('0   ' + str(spinMult) + '\n')
        for idx in range(0, len(buildTS.atoms)):
            atSym = buildTS.rdMol.GetAtomWithIdx(idx).GetSymbol()
            position = buildTS.rdMol.GetConformer(buildTS.rdMolConfId).GetAtomPosition(idx)
            xPt = "  " + "%.7f" %position.x
            yPt = "  " + "%.7f" %position.y
            zPt = "  " + "%.7f" %position.z
            if position.x >= 0.0:
                xPt = " " + xPt
            if position.y >= 0.0:
                yPt = " " + yPt
            if position.z >= 0.0:
                zPt = " " + zPt
            fout.write(atSym + xPt + yPt + zPt + "\n")
        fout.write("\n")
        fout.close()
        
        fout = open(os.path.join(newPath, title+'.sh'), "w")
        fout.write('#!/bin/sh\n')
        fout.write('#BSUB -q normal\n')
        fout.write('#BSUB -o ./output/' + title + '.out\n')
        fout.write('#BSUB -J reactionTS\n\n')
        fout.write('export GAUSS_EXEDIR=/share/apps/g09\n')
        fout.write('export PATH=$GAUSS_EXEDIR:$PATH\n\n')
        fout.write('g09 < ' + filename + ' > ./log/' + title + '.log\n\n')
        fout.close()