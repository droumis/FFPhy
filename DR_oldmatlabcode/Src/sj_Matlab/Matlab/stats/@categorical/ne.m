%NE Not equal for categorical arrays.
%   TF = NE(A,B) returns a logical array the same size as the categorical
%   arrays A and B, containing true (1) where the corresponding elements of A
%   and B are not equal, and false (0) otherwise.  A and B need not have the
%   same sets of levels.  The test is performed by comparing the labels of the
%   levels of each pair of elements.
%
%   TF = NE(A,S) or TF = NE(S,B), where S is a character string, returns a
%   logical array the same size as A or B, containing true (1) where the
%   levels of corresponding elements of A or B have labels not equal to S.
%
%   Elements with undefined levels are considered not equal to any other value,
%   including other undefined elements.
%
%   See also CATEGORICAL/EQ.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/12/15 19:30:57 $

% This is an abstract method.