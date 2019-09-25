function r=gammar(m,n,a);
%GAMMAR Generates gamma random deviates.
%	GAMMAR(M,N,A) is an MxN matrix of random deviates from the standard
%	gamma distribution with shape parameter A.  A must be a scalar.
%
%	B*GAMMAR(M,N,A) is a matrix of random deviates from the gamma
%	distribution with shape parameter A and scale parameter B.  The
%	distribution then has mean A*B and variance A*B^2.
%
%	See GAMMAP, GAMMAQ, RAND.

%  Gordon Smyth, University of Queensland, gks@maths.uq.edu.au
%  9 Dec 1999.

r = zeros(m,n);
for i=1:m,
   for j=1:n,
      r(i,j)=gammar1(a);
   end
end
