function a = reorderlevels(a,newlevels)
%REORDERLEVELS Reorder levels in a categorical array.
%   B = REORDERLEVELS(A,NEWLEVELS) reorders the levels of the categorical array A.
%   NEWLEVELS is a cell array of strings or a 2-dimensional character matrix
%   that specifies the new order.  NEWLEVELS must be a reordering of GETLABELS(A).
%
%   The order of the levels of a categorical array has no significance other
%   than for display purposes.
%
%   See also CATEGORICAL/ADDLEVELS, CATEGORICAL/DROPLEVELS, CATEGORICAL/ISLEVEL,
%            CATEGORICAL/MERGELEVELS, CATEGORICAL/GETLABELS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:00 $

newlevels = checklabels(newlevels,2); % error if duplicates

[tf,convert] = ismember(a.labels,newlevels);
if (length(newlevels) ~= length(a.labels)) || ~all(tf)
    error('stats:categorical:reorderlevels:InvalidNewLevels', ...
          'NEWLEVELS must be a permutation of the existing categorical levels.');
end
convert = cast([0 convert],class(a.codes)); % there may be zeros in a.codes
a.codes = reshape(convert(a.codes+1), size(a.codes));
a.labels = newlevels;
