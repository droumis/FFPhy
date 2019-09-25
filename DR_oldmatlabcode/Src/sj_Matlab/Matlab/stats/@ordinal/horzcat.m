function a = horzcat(varargin)
%HORZCAT Horizontal concatenation for ordinal arrays.
%   C = HORZCAT(A, B, ...) horizontally concatenates the ordinal arrays A, B,
%   ... .  For matrices, all inputs must have the same number of rows.  For
%   N-D arrays, all inputs must have the same sizes except in the second
%   dimension.  All inputs must have the same sets of ordinal levels,
%   including their order.
%
%   C = HORZCAT(A,B) is called for the syntax [A B].
%
%   See also ORDINAL/CAT, ORDINAL/VERTCAT.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:28 $

a = varargin{1};

for i = 2:nargin
    b = varargin{i};
    if ~isa(b,class(a))
        error('stats:ordinal:horzcat:TypeMismatch', ...
              'All input arguments must be from the same categorical class.');
    end
    if isequal(a.labels,b.labels)
        bcodes = b.codes;
    else
        error('stats:ordinal:horzcat:OrdinalLevelsMismatch', ...
              'Ordinal levels and their ordering must be identical.');
    end

    try
        a.codes = horzcat(a.codes, bcodes);
    catch
        rethrow(lasterror);
    end
end
