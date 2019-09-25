function a = vertcat(varargin)
%VERTCAT Vertical concatenation for nominal arrays.
%   C = VERTCAT(A, B, ...) vertically concatenates the nominal arrays A,
%   B, ... .  For matrices, all inputs must have the same number of columns.
%   For N-D arrays, all inputs must have the same sizes except in the first
%   dimension.  The set of nominal levels for C is the sorted union of the
%   sets of levels of the inputs, as determined by their labels.
%
%   C = VERTCAT(A,B) is called for the syntax [A; B].
%
%   See also NOMINAL/CAT, NOMINAL/HORZCAT.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:20 $

a = varargin{1};

for i = 2:nargin
    b = varargin{i};
    if ~isa(b,class(a))
        error('stats:nominal:vertcat:TypeMismatch', ...
              'All input arguments must be from the same categorical class.');
    end
    if isequal(a.labels,b.labels)
        bcodes = b.codes;
    else
        % Get a's codes for b's data, possibly adding to a's levels
        [bcodes,a.labels] = matchlevels(a,b);
    end

    try
        a.codes = vertcat(a.codes, bcodes);
    catch
        rethrow(lasterror);
    end
end
