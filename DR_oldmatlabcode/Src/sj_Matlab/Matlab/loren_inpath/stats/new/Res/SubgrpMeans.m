% SubgrpMeans: Given two vectors of group identifiers for groups and subgroups (as would be
%         used in a nested anova or analysis of homogeneous subsets), produces a table of
%         sample sizes and cell means (for a single variable).
%
%     Usage: [N,cellmeans,ugrp] = subgrpmeans(x,grps,subgrps)
%
%         x =         data vector.
%         grps =      group-identification vector for k groups.
%         subgrps =   subgroup-identification vector for max m subgroups/group.
%         ---------------------------------------------------------------------
%         N =         [k x m] matrix of sample sizes per cell.
%         cellmeans = [k x m] matrix of corresponding cell means.
%         ugrp =      [k x 1] vector of sorted group identifiers, corresponding
%                       to rows of output matrices.
%

% RE Strauss, 2/21/02
%   2/26/02 - initialize N before loop.

function [N,cellmeans,ugrp] = subgrpmeans(x,grps,subgrps)
  if (~isvector(x))
    error('  SubgrpMeans: data matrix must be a single vector.');
  end;
  
  ugrp = uniquef(grps,1);
  N = [];
  for ig = 1:length(ugrp)
    i = find(grps == ugrp(ig));
    sg = subgrps(i);
    xg = x(i);
    [usubgrp,fsubgrp] = uniquef(sg);
    fsubgrp = fsubgrp';
    [N,fsubgrp] = padcols(N,fsubgrp);
    N(ig,:) = fsubgrp;
    for isg = 1:length(usubgrp)
      j = find(sg == usubgrp(isg));
      cellmeans(ig,isg) = mean(xg(j));
    end;
  end;

  return;
  
