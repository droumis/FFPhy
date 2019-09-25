function [out] = JY_getreroutetrialvelocity(index, excludetimes, occupancyoutput)
% gets the linear velocity


warning('OFF','MATLAB:divideByZero');

%linpos = linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells;
%data = data{index(1)}{index(2)}.Pos.correcteddata
if ~isempty(occupancyoutput{index(1)}{index(2)})
    
    
    out.occupancy=occupancyoutput{index(1)}{index(2)}.occupancy;
    
else
   
     out.occupancy=[];
end


warning('ON','MATLAB:divideByZero');