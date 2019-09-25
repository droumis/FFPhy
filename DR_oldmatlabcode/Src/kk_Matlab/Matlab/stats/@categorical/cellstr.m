function b = cellstr(a)
%CELLSTR Convert categorical array to cell array of strings.
%   B = CELLSTR(A) converts the categorical array A to a cell array of
%   strings.  Each element of B contains the categorical level label for the
%   corresponding element of A.
%
%   See also CATEGORICAL/CHAR, CATEGORICAL/GETLABELS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:25 $

labs = [categorical.undefLabel a.labels];
b = reshape(labs(a.codes+1),size(a.codes));
