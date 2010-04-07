function [dataArray_RMG] = convertNASA2RMG(dataArray_polys)
%convertNASA2RMG.m
    
[rows,columns] = size(dataArray_polys);
for n=1:rows
    T_lower(1) = dataArray_polys{n,2};
    T_upper(1) = dataArray_polys{n,4};
    H(1,:) = [dataArray_polys{n,6} dataArray_polys{n,7}/2 ...
        dataArray_polys{n,8}/3 dataArray_polys{n,9}/4 ...
        dataArray_polys{n,10}/5 dataArray_polys{n,11}];
    S(1,:) = [dataArray_polys{n,6} dataArray_polys{n,7} ...
        dataArray_polys{n,8}/2 dataArray_polys{n,9}/3 ...
        dataArray_polys{n,10}/4 dataArray_polys{n,12}];
    Cp(1,:) = [dataArray_polys{n,6} dataArray_polys{n,7} ...
        dataArray_polys{n,8} dataArray_polys{n,9} ...
        dataArray_polys{n,10}];
    
    T_lower(2) = dataArray_polys{n,13};
    T_upper(2) = dataArray_polys{n,15};
    H(2,:) = [dataArray_polys{n,17} dataArray_polys{n,18}/2 ...
        dataArray_polys{n,19}/3 dataArray_polys{n,20}/4 ...
        dataArray_polys{n,21}/5 dataArray_polys{n,22}];
    S(2,:) = [dataArray_polys{n,17} dataArray_polys{n,18} ...
        dataArray_polys{n,19}/2 dataArray_polys{n,20}/3 ...
        dataArray_polys{n,21}/4 dataArray_polys{n,23}];
    Cp(2,:) = [dataArray_polys{n,17} dataArray_polys{n,18} ...
        dataArray_polys{n,19} dataArray_polys{n,20} ...
        dataArray_polys{n,21}];
    
    dataArray_RMG{n,1} = dataArray_polys{n,1};
    dataArray_RMG{n,2} = dataArray_polys{n,24};
    foundValidTRange = false;

    for m=1:2
        if T_lower(m) <= 298 && 298 <= T_upper(m)
            dataArray_RMG{n,3} = H(m,:) * [1 298 298^2 298^3 298^4 1/298]' * ...
                1.987e-3 * 298;
            dataArray_RMG{n,4} = S(m,:) * [log(298) 298 298^2 298^3 298^4 1]' ...
                * 1.987;
            foundValidTRange = true;
            break;
        end
    end
    
    if ~foundValidTRange
        diff = [T_lower(1)-298 T_lower(2)-298];
        index = find(diff==min(diff));
        dataArray_RMG{n,3} = H(index,:) * [1 298 298^2 298^3 298^4 1/298]' * ...
                1.987e-3 * 298;
        dataArray_RMG{n,4} = S(index,:) * [log(298) 298 298^2 298^3 298^4 1]' ...
                * 1.987;
        dataArray_RMG{n,12} = 'Neither TRange included 298K!';
    end

    Temps = [300 400 500 600 800 1000 1500];
    for numCp=1:7
        T = [1 Temps(numCp) Temps(numCp)^2 Temps(numCp)^3 Temps(numCp)^4]';
        for m=1:2
            if T_lower(m) <= Temps(numCp) && Temps(numCp) < T_upper(m)
                dataArray_RMG{n,numCp+4} = Cp(m,:) * T * 1.987;
                break;
            end
        end
    end
end
    
return