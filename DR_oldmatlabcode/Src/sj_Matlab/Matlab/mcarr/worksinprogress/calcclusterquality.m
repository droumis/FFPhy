function out = calcclusterquality(index, excludetimes, clust)
% Extracts the cluster quality from all cells listed in index
% Out is an N x 2 vector
%   N = number of cells
%       1st column: isolation distance
%       2nd column: lratio


out = [clust{index(1)}{index(2)}{index(3)}{index(4)}.isoldist clust{index(1)}{index(2)}{index(3)}{index(4)}.lratio];

end
