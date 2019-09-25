function t = isequal(a,varargin)
%ISEQUAL True if categorical arrays are equal.
%   TF = ISEQUAL(A,B) is true (1) if the categorical arrays A and B are the
%   same class, have the same size and the same sets of levels, and contain
%   the same values, and false (0) otherwise.
%
%   TF = ISEQUAL(A,B,C,...) is true (1) if all the input arguments are equal.
%
%   Elements with undefined level are not considered equal to each other.
%
%   See also CATEGORICAL/GETLABELS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:47 $


if nargin < 2
    error('Not enough input arguments.');
end
t = isa(a,'categorical') && all(a.codes(:) > 0);
for i = 1:length(varargin)
    if ~t, break; end
    b = varargin{i};
    if isa(b,class(a)) && all(b.codes(:) > 0)
        t = isequal(a.codes,b.codes) && isequal(a.labels,b.labels);
    else
        t = false;
    end
end
