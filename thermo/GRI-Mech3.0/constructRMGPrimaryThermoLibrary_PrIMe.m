function constructRMGPrimaryThermoLibrary_PrIMe
%constructRMGPrimaryThermoLibrary_PrIMe.m
%   Function generates a PrimaryThermoLibrary for RMG software.  The
%       two files generated are primeID.txt (a list of primeIDs and their
%       respective InChIs) and Library.txt.  The primeID.txt file may be
%       fed to the RMG module InChI2AdjList to convert the InChIs to
%       adjacency lists.  The output file from this module is the
%       Dictionary.txt file required for the primaryThermoLibrary folder.
%
% Written by MRH on 30Jun2009

% Search the PrIMe species catalog for recommended thermochemistry data
recommendedThermo = extractRecommendedPolys;
% Convert the NASA-7 polynomials to RMG thermo format (Hf298, S298, etc.)
rmgThermo = convertNASA2RMG(recommendedThermo);
% Write the primeID.txt and Library.txt files
writePrimaryThermoLibraryFiles(rmgThermo);

return