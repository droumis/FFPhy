function b = char(a)
%CHAR Convert categorical array to character array.
%   B = CHAR(A) converts the categorical array A to a 2-dimensional character
%   matrix.  CHAR does not preserve the shape of A.  B contains NUMEL(A) rows,
%   and each row of B contains the categorical level label for the
%   corresponding element of A(:).
%
%   See also CATEGORICAL/CELLSTR, CATEGORICAL/GETLABELS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:26 $

labs = [categorical.undefLabel a.labels];
b = char(labs(a.codes+1));
