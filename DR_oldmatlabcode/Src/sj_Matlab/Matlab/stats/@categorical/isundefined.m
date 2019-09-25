function tf = isundefined(a)
%ISUNDEFINED True for elements of a categorical array that are undefined.
%   TF = ISUNDEFINED(A) returns a logical array the same size as the
%   categorical array A, containing true (1) where the corresponding element
%   of A is undefined, i.e., does not have a value from the A's set of
%   categorical levels, and false (0) otherwise.
%
%   See also CATEGORICAL/ISMEMBER.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:51 $

tf = (a.codes == 0);
