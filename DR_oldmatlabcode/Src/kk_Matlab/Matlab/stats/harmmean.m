function m = harmmean(x,dim)
%HARMMEAN Harmonic mean.
%   M = HARMMEAN(X) returns the harmonic mean of the values in X.  For
%   vector input, M is the inverse of the mean of the inverses of the
%   elements in X.  For matrix input, M is a row vector containing the
%   harmonic mean of each column of X.  For N-D arrays, HARMMEAN operates
%   along the first non-singleton dimension.
%
%   HARMMEAN(X,DIM) takes the harmonic mean along dimension DIM of X.
%
%   See also MEAN, GEOMEAN, TRIMMEAN.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.8.2.3 $  $Date: 2004/07/28 04:38:28 $

% Take the reciprocal of the mean of the reciprocals of X.
if nargin == 1
    m = 1 ./ mean(1./x);
else
    m = 1 ./ mean(1./x, dim);
end
