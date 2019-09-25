function [resweepn, map] = collapse_sweeps(sweepn)

% COLLAPSE SWEEPS: Renumbers sweep assignments for rasters: To get a continuous Y-axis


sweeps = unique(sweepn);  % get a list of unique labels . . .
numsweeps = length(sweeps);    %

resweepn = zeros(size(sweepn));
for n = 1:numsweeps
    resweepn(find(sweepn == sweeps(n,1))) = n;
end
map = [unique(sweepn) unique(resweepn)]; 

return;
