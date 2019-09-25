function v=cutpoint(t,j)
%CUTPOINT Cutpoints used for branches in decision tree.
%   V=CUTPOINT(T) returns an N-element vector V of the values used as
%   cutpoints in the decision tree T.  For each branch node J based on a
%   continuous predictor variable Z, the left child is chosen if Z<V(J) and
%   the right child is chosen if Z>=V(J).  V is NaN for branch nodes based
%   on categorical predictors and for non-branch (leaf) nodes.
%
%   V=CUTPOINT(T,J) takes an array J of node numbers and returns the
%   cutpoints for the specified nodes.
%
%   See also CLASSREGTREE, CLASSREGTREE/CUTVAR, CLASSREGTREE/CUTCATEGORIES, CLASSREGTREE/CUTTYPE.

%   Copyright 2006-2007 The MathWorks, Inc. 
%   $Revision: 1.1.6.1.2.2 $  $Date: 2007/01/30 02:22:30 $

if nargin>=2 && ~validatenodes(t,j)
    error('stats:classregtree:parent:InvalidNode',...
          'J must be an array of node numbers or a logical array of the proper size.');
end

% Get variable numbers and cut points
if nargin<2
    n = t.var;
    v = t.cut;
else
    n = t.var(j,:);
    v = t.cut(j,:);
end

% Remove invalid values
v(n<=0) = NaN;