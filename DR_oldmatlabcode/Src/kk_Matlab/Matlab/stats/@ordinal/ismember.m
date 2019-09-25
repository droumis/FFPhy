function [tf,loc] = ismember(a,s)
%ISMEMBER True for elements of an ordinal array in a set.
%   TF = ISMEMBER(A,LEVELS) returns a logical array the same size as the
%   ordinal array A, containing true (1) where the level of the corresponding
%   element of A is equal to one of the levels specified in LEVELS, and false
%   (0) otherwise.  LEVELS is an ordinal array, or a cell array of strings or
%   2-dimensional character array containing level labels.
%
%   [TF,LOC] = ISMEMBER(A,LEVELS) also returns an index array LOC containing
%   the highest absolute index in LEVELS for each element in A whose level is
%   a member of LEVELS, and 0 if there is no such index.
%
%   See also ORDINAL/ISLEVEL, ORDINAL/UNIQUE, ORDINAL/UNION,
%            ORDINAL/INTERSECT, ORDINAL/SETXOR, ORDINAL/SETDIFF.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:32:30 $

if isa(s,'categorical')
    if ~isa(s,class(a))
        error('stats:ordinal:ismember:TypeMismatch', ...
              'All input arguments must be from the same categorical class.');
    elseif ~isequal(a.labels,s.labels)
        error('stats:ordinal:iismember:OrdinalLevelsMismatch', ...
              'Ordinal levels and their ordering must be identical.');
    end
    s = cellstr(s(:));
end

[dum,scodes] = ismember(s,a.labels);

if nargout < 2
    tf = ismember(a.codes,scodes(scodes>0));
else
    [tf,loc] = ismember(a.codes,scodes(scodes>0));
end