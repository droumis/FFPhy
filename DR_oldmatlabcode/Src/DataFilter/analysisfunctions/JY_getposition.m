function [out] = JY_getreroutetrialvelocity(index, excludetimes, data, linpos)
% the positions on the track for all excluded periods


warning('OFF','MATLAB:divideByZero');

%linpos = linpos{index(1)}{index(2)}.statematrix.linearDistanceToWells;
%data = data{index(1)}{index(2)}.Pos.correcteddata
if ~isempty(linpos{index(1)}{index(2)})
    

  % get all positions during ripples ie. the excluded times
  out.rippleposition=data{index(1)}{index(2)}.Pos.correcteddata(~isExcluded(data{index(1)}{index(2)}.Pos.correcteddata(:,1), excludetimes),2:3);
  
  % get all velocities during ripples
  
  out.ripplevelocity=linpos{index(1)}{index(2)}.statematrix.linearVelocity(~isExcluded(data{index(1)}{index(2)}.Pos.correcteddata(:,1), excludetimes),1);
  
  
else

    out.rippleposition=[];
    out.ripplevelocity=[];
end

out.allpositions=data{index(1)}{index(2)}.Pos.correcteddata;
warning('ON','MATLAB:divideByZero');