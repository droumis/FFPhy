function out=diagsum(X)
%DIAGSUM Sum of diagonals of a matrix.
%	DIAGSUM(A) is a row vector containing the sums of the diagonals
%	of the matrix A, from lower left to upper right.

% From makra@athena.mit.edu (Mohamad A. Akra)
% Massachusetts Institute of Technology
% Posted to comp.soft-sys.matlab, 29 July 93

[M,N] = size(X);
if M == 1,
   out = X;
else
   Y = [X zeros(M,M-2)].';
   out = sum( reshape( [zeros(M-1,1); Y(:); 0], M+N-1, M ).');
end;
