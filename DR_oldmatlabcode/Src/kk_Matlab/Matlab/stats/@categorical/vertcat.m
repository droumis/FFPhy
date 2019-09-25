%VERTCAT Vertical concatenation for categorical arrays.
%   C = VERTCAT(A, B, ...) vertically concatenates the categorical arrays A,
%   B, ... .  For matrices, all inputs must have the same number of columns.
%   For N-D arrays, all inputs must have the same sizes except in the first
%   dimension.  The set of categorical levels for C is the sorted union of the
%   sets of levels of the inputs, as determined by their labels.
%
%   C = VERTCAT(A,B) is called for the syntax [A; B].
%
%   See also CATEGORICAL/CAT, CATEGORICAL/HORZCAT.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/12/15 19:31:23 $

% This is an abstract method.