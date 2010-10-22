function [GoldStarThermo] = extractRecommendedPolys
%writePrimaryThermoLibrary.m
%   Function searches the PrIMe species catalog for species with PrIMe-
%       recommended thermodynamic polynomials (thp00000000.xml).  For those
%       species that have recommended data, the data file (thp0000000?.xml)
%       is opened and the polynomials (with their lower and upper
%       temperature bounds) are extracted and stored in the cell array
%       GoldStarThermo.
%
%   NOTE: The ReactionLab GUI must be opened for this function to run
%       >> primekinetics
%
%   Written by MRH on 30-Jun-2009

numGSSpecies = 0;
for counter = 0:10889
    foundRecThermo = true;
    %Attempt to grab the link to the PrIMe recommended thermo file
    try
        primeThermoObj0 = gate2primeData('get',...
            {int2PrIMeID(counter,'s'),'thp00000000'});
    catch
        foundRecThermo = false;
    end
    if (foundRecThermo && ~isempty(primeThermoObj0))
        numGSSpecies = numGSSpecies + 1;
        %Determine which thermo data file is the recommended file
        recThermo = getField(primeThermoObj0.node, ...
            'speciesLink.thermodynamicPolynomialsLink._primeID');
        %Grab the actual recommended thermo file
        primeThermoObj = gate2primeData('get',...
            {int2PrIMeID(counter,'s'),char(recThermo)});
        %Grab all thermo polynomials and their valid range
        thermoPolys = getPolynomials(primeThermoObj);
        %Fill a cell array with thermodynamic information for each species
        GoldStarThermo{numGSSpecies,1} = int2PrIMeID(counter,'s');
        %Loop over the number of thermo polynomial sets
        for numPolySets = 1:length(thermoPolys)
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+2} = ...
                thermoPolys(numPolySets).Tmin.value;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+3} = ...
                thermoPolys(numPolySets).Tmin.units;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+4} = ...
                thermoPolys(numPolySets).Tmax.value;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+5} = ...
                thermoPolys(numPolySets).Tmax.units;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+6} = ...
                thermoPolys(numPolySets).poly.a1;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+7} = ...
                thermoPolys(numPolySets).poly.a2;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+8} = ...
                thermoPolys(numPolySets).poly.a3;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+9} = ...
                thermoPolys(numPolySets).poly.a4;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+10} = ...
                thermoPolys(numPolySets).poly.a5;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+11} = ...
                thermoPolys(numPolySets).poly.a6;
            GoldStarThermo{numGSSpecies,(numPolySets-1)*11+12} = ...
                thermoPolys(numPolySets).poly.a7;
        end
        %Grab the PrIMe species obj
        primeSpeciesObj = gate2primeData('get',{int2PrIMeID(counter,'s')});
        GoldStarThermo{numGSSpecies,24] = getInChI(primeSpeciesObj);
    end
end

return

function [PrIMeID] = int2PrIMeID(integer, primeIDtype)
%int2PrIMeID.m

PrIMeID = primeIDtype;
if (integer >= 0 && integer < 1e8)
    if (integer < 1e1)
        PrIMeID = strcat(PrIMeID, '0000000', int2str(integer));
    elseif (integer < 1e2)
        PrIMeID = strcat(PrIMeID, '000000', int2str(integer));
    elseif (integer < 1e3)
        PrIMeID = strcat(PrIMeID, '00000', int2str(integer));
    elseif (integer < 1e4)
        PrIMeID = strcat(PrIMeID, '0000', int2str(integer));
    elseif (integer < 1e5)
        PrIMeID = strcat(PrIMeID, '000', int2str(integer));
    elseif (integer < 1e6)
        PrIMeID = strcat(PrIMeID, '00', int2str(integer));
    elseif (integer < 1e7)
        PrIMeID = strcat(PrIMeID, '0', int2str(integer));
    else
        PrIMeID = strcat(PrIMeID, int2str(integer));
    end
else
    disp('Error: value2PrIMeID.m requires an integer >= 0 and < 1e8');
end

return;