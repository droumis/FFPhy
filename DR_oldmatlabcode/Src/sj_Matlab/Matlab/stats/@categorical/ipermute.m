function b = ipermute(a,order)
%IPERMUTE Inverse permute dimensions of a categorical array.
%   A = IPERMUTE(B,ORDER) is the inverse of PERMUTE. IPERMUTE rearranges the
%   dimensions of the categorical array B so that PERMUTE(A,ORDER) will
%   produce B.  The array produced has the same values of A but the order of
%   the subscripts needed to access any particular element are rearranged as
%   specified by ORDER.  The elements of ORDER must be a rearrangement of the
%   numbers from 1 to N.
%
%   See also CATEGORICAL/PERMUTE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:45 $

b = a;
b.codes = ipermute(a.codes,order);
