function a = mat2cell(m, rows, cols)
% mat2cell - burst matrix into cell array of submatrices
%
%	A = MAT2CELL(M, R, C) bursts the matrix M into submatrices
%	described by R and C and returns a cell array of the
%	corresponding matrices.  R and C may be vectors, containing
%	the row or column indices that begin the submatrices.
%	Alternatively, they may be scalars, in which case the rows or
%	columns of M are divided equally among the submatrices.
%
%	See also NUM2CELL

% Copyright (c) 1997 Maneesh Sahani (maneesh@caltech.edu)

error(nargchk(3,3,nargin));

[r, c] = size(m);

if length(rows) == 1
  rows = 1:floor(r/rows):r;
end

if length(cols) == 1
  cols = 1:floor(c/cols):c;
end

a = cell(length(rows), length(cols));

rows = [rows, r+1];
cols = [cols, c+1];

for i = 1:length(rows)-1
  for j = 1:length(cols)-1
    a{i,j} = m(rows(i):rows(i+1)-1, cols(j):cols(j+1)-1);
  end
end  
