function element = index(matrix, i, j)
% index - select elements from a matrix
%
%	SUB = INDEX(MAT, I, J) returns the submatrix defined by
%	MAT(I,J).  This is useful to extract submatrices from the
%	return values of function, since matlab does not allow the
%	construct foo(bar)(1,1).

% Copyright 1997 John Pezaris (pz@caltech.edu), 
%		 Maneesh Sahani (maneesh@caltech.edu) 

if (nargin < 3)
   element = matrix(i);
else
   element = matrix(i, j);
end

