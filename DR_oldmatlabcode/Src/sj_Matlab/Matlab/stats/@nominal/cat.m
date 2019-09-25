function a = cat(dim,varargin)
%CAT Concatenate nominal arrays.
%   C = CAT(DIM, A, B, ...) concatenates the nominal arrays A, B, ...
%   along dimension DIM.  All inputs must have the same size except along
%   dimension DIM.  The set of nominal levels for C is the sorted union of
%   the sets of levels of the inputs, as determined by their labels.
%
%   See also NOMINAL/HORZCAT, NOMINAL/VERTCAT.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:09 $

a = varargin{1};

for i = 2:nargin-1
    b = varargin{i};
    if ~isa(b,class(a))
        error('stats:nominal:cat:TypeMismatch', ...
              'All input arguments must be from the same categorical class.');
    end
    if isequal(a.labels,b.labels)
        bcodes = b.codes;
    else
        % Get a's codes for b's data, possibly adding to a's levels
        [bcodes,a.labels] = matchlevels(a,b);
    end

    try
        a.codes = cat(dim, a.codes, bcodes);
    catch
        rethrow(lasterror);
    end
end
