function [c,varargout] = max(a,b,dim)
%MAX Largest element in an ordinal array.
%   B = MAX(A), when A is an ordinal vector, returns the largest element in A.
%   For ordinal matrices, MAX(A) is a row vector containing the maximum
%   element from each column.  For N-D ordinal arrays, MAX(A) operates along
%   the first non-singleton dimension.  B is an ordinal array with the same
%   levels as A.
%
%   [B,I] = MAX(A) returns the indices of the maximum values in vector I.
%   If the values along the first non-singleton dimension contain more
%   than one maximal element, the index of the first one is returned.
%
%   C = MAX(A,B) returns an ordinal array the same size as A and B with the
%   largest elements taken from A or B.  A and B must have the same sets of
%   ordinal levels, including their order.
%
%   [B,I] = MAX(A,[],DIM) operates along the dimension DIM. 
%
%   Elements with undefined levels are not considered greater than any other
%   elements, including each other.
%
%   See also ORDINAL/MIN, ORDINAL/SORT.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:34 $

if nargin > 2 && ~isnumeric(dim)
    error('stats:ordinal:max:InvalidDim', ...
          'DIM must be a positive integer scalar in the range 1 to 2^31.');
end

if (nargin < 2) || isequal(b,[]) % && isa(a,'ordinal')
    acodes = a.codes;
    bcodes = [];
    c = ordinal([],a.labels,1:length(a.labels));
else
    [acodes,bcodes] = ordinalcheck(a,b);
    if ischar(a)
        c = ordinal([],b.labels,1:length(b.labels));
    else
        c = ordinal([],a.labels,1:length(a.labels));
    end
end

% Undefined elements have code zero, less than any legal code.  They will not
% be the max value unless there's nothing else.
try
    if nargin == 1
        [ccodes,varargout{1:nargout-1}] = max(acodes);
    elseif nargin == 2
        [ccodes,varargout{1:nargout-1}] = max(acodes,bcodes);
    else
        [ccodes,varargout{1:nargout-1}] = max(acodes,bcodes,dim);
    end
catch
    rethrow(lasterror);
end
c.codes = ccodes;
