function writePrimaryThermoLibraryFiles(dataArray_RMG)
%writePrimaryThermoLibraryFiles.m

%Open a text file
fidLibrary = fopen('Library.txt','wt');
fidInChI = fopen('InChI.txt','wt');

fprintf(fidLibrary,'//Library.txt\n' + ...
    '//\tversion 1.0\n\n' + ...
    '//This file contains PrIMe-recommended thermochemistry data ' + ...
    'as of ' + datestr(now) + '\n' + ...
    '//All data was converted from NASA-7 polynomial format' + ...
    ' (unless otherwise noted)\n' + ...
    '//All uncertainty values were assigned values of 0.0\n\n' + ...
    '//Name\tHf298\tS298\tCp300\tCp400\tCp500\tCp600\tCp800\t' + ...
    'Cp1000\tCp1500\tdH\tdS\tdCp\tNotes\n\n' + ...
    '//Units\n' + ...
    '//H: kcal/mol\n' + ...
    '//S and Cp: cal/mol/K\n\n');
%Loop over all "Gold Star" species
sizeCell = size(dataArray_RMG);
for numGSSpecies = 1:sizeCell(1)
    %The primeID string (for now) will be the "Group" name
    fprintf(fidLibrary,dataArray_RMG{numGSSpecies,1} + '\t');
    %Fill in the thermodynamic (H, S, Cp's) data
    for numParams = 1:9
        fprintf(fidLibrary,...
            num2str(dataArray_RMG{numGSSpecies,numParams+2}) + '\t');
    end
    %Fill in the dH, dS, and dCP values (set to zero)
    fprintf(fidLibrary,'0.0\t0.0\t0.0\t');
    %Comment where the data came from
    if ~isempty(dataArray_RMG{numGSSpecies,12})
        fprintf(fidLibrary,dataArray_RMG{numGSSpecies,12});
    end
    fprintf(fidLibrary,'\n');
    
    fprintf(fidInChI,dataArray_RMG{numGSSpecies,1} + '\t' + ...
        dataArray_RMG{numGSSpecies,2} + '\n');
end

fclose(fidLibrary);
fclose(fidInChI);

return