
function out = DFAsj_getclustqual(index, excludetimes, clustqual, varargin)

if ~isempty(excludetimes)
    excludetimes = [];
end

% Put in structure and return
% -----------------------------

out.index = index;
%  if index(1)==5 & index(2)==1 & index(3)==4
%      keyboard;
%  end
% index
if ~isempty(clustqual{index(1)})
    
    %if ~isempty(clustqual{index(1)}{index(2)}{index(3)}{index(4)})
    out.isoldist = clustqual{index(1)}{index(2)}{index(3)}{index(4)}.isoldist;
    out.lratio = clustqual{index(1)}{index(2)}{index(3)}{index(4)}.lratio;
    %end
else
    out.isoldist = NaN;
    out.lratio = NaN;
end
