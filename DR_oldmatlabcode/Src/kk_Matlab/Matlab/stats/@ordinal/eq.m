function t = eq(a,b)
%EQ Equality for ordinal arrays.
%   TF = EQ(A,B) returns a logical array the same size as the ordinal arrays A
%   and B, containing true (1) where the corresponding elements of A and B are
%   equal, and false (0) otherwise.  A and B must have the same sets of
%   ordinal levels, including their order.
%
%   TF = EQ(A,S) or TF = EQ(S,B), where S is a character string, returns a
%   logical array the same size as A or B, containing true (1) where the
%   elements of A or B have levels whose labels are equal to S.
%
%   Elements with undefined levels are not considered equal to each other.
%
%   See also ORDINAL/NE, ORDINAL/GE, ORDINAL/LE, ORDINAL/LT,  ORDINAL/GT.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:24 $

[acodes,bcodes] = ordinalcheck(a,b);

% undefined elements cannot be equal to anything
t = (acodes == bcodes) & (acodes ~= 0);
