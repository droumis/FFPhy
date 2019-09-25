function f = levels(nlevels,nreps,length)
%LEVELS    Generate factor levels.  LEVELS(nlevels,nreps,length) is
%          kron( 1:nlevels, ones(nreps,1) ) repeated to make a vector
%          of length elements.

n=ceil( length/(nlevels*nreps) );
f=kron( ones(n,1) , kron( (1:nlevels)' , ones(nreps,1) ) );
f=f(1:length);
