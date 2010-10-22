7-Jul-2009 (MRH)

folder: $RMG/databases/RMG_database/thermo/GRI-Mech3.0

This folder contains the following files (in addition to this readme.txt file):

- constructRMGPrimaryThermoLibrary_PrIMe.m
	This is the main MATLAB function.  Running this program will mine the PrIMe Warehouse species
		catalog for all recommended thermodynamic polynomials (thp00000000.xml) files.  Note that
		the PrIMe GUI ReactionLab must be open for this program to run properly (invoked by typing
		>> primekinetics at the Command Window).

- extractRecommendedPolys.m
	This function is called by the main MATLAB function.  The purpose of this function is to mine
		the PrIMe species catalog for all species with thp00000000.xml files and extract the 
		important information.  This function requires no inputs and outputs a cell array, 
		containing the following information for each species in the PrIMe species catalog that 
		has recommended thermodynamic polynomials:
	* PrIMe species ID
	* NASA-7 formatted thermodynamic polynomials, including the lower and upper temperature bound
		(w/units)
	* InChI string (if available)
	An example output from this function is located in the cellArrays.mat file: GoldStarThermo

- convertNASA2RMG.m
	This function is called by the main MATLAB function.  The purpose of this function is to 
		convert the NASA-7 polynomials into RMG-readable thermodynamic data.  This function 
		requires the cell array	output from extractRecommendedPolys.m as input and outputs its 
		own cell array, containing the following information for each each species in the PrIMe 
		species catalog that has recommended thermodynamic polynomials:
	* PrIMe species ID
	* InChI string (if available)
	* Hf298 (kcal/mol)
	* S298 (cal/mol/K)
	* Cp300, Cp400, Cp500, Cp600, Cp800, Cp1000, Cp1500 (cal/mol/K)
	* Other information (e.g. When neither of the NASA-7 polynomial's valid temperature range
		contains 298K, troublesome for computing Hf298 and S298)
	An example output from this function is located in the cellArray.mat file: GoldStarRMG

- writePrimaryThermoLibraryFiles.m
	This function is called by the main MATLAB function.  The purpose of this function is to
		write the RMG-required Library.txt and Dictionary.txt files that must be present for
		any PrimaryThermoLibrary.  This function requires the cell array output from 
		convertNASA2RMG.m as input; there is no output, but two files (Library.txt and InChI.txt)
		are written.

- cellArrays.mat
	File contains 2 cell arrays, GoldStarThermo and GoldStarRMG.  These are example outputs from
		the functions extractRecommendedPolys.m and convertNASA2RMG.m, respectively.

- Library.txt
	This file contains a list of species and their RMG-readable thermodynamic values.  A species
		is labeled with its PrIMe species ID (e.g. s00009015), as this is PrIMe's unique
		species identifier.

- InChI.txt
	This file contains a list of species and their InChI strings.  A species is labeled with its
		PrIMe species ID.  This file is fed to the RMG module InChI2AdjList.java, present in the
		default package of the RMG software, in order to construct an adjacency list for each
		species present in Library.txt

- adjList_output.txt
	This file is the output from feeding InChI.txt to InChI2AdjList.java.  It contains a list of
		species, identified by their PrIMe species ID, and their adjacency list.  This is the
		exact output from the RMG module.

- Dictionary.txt
	This file contains the same information as adjList_output.txt, with some minor corrections
		by MRH.  The corrections include:
	* Commenting out all Nitrogen-containing species (since RMG does not recognize N)
	* Commenting out the C3H7 species (since the GRI-Mech 3.0 C3H7 species is a lumped species
		of 1-C3H7 and 2-C3H7)
	* Commenting out Argon (since RMG does not recognize Ar)
	* Correcting some of the adjacency lists, e.g. CH3, CH2(S), CH2, OH, and O.  The InChI
		executable does not label these as radicals (using the M  RAD tag near the bottom of
		the .mol file).  Rather, it labels the bonds to the heavy atoms.  MRH has not updated
		the Java source code to handle this yet.

Together, Dictionary.txt and Library.txt make up a new PrimaryThermoLibrary: GRI-Mech3.0