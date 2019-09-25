function ind = indicato(x,setto0)
%INDICATO INDICATO(X) constructs indicator variables for the levels
%	taken by the vector or matrix X.  Levels are rounded down to
%	integers.  INDICATO(X,SETTO0) excludes the SETTO0's column.

%	GKS  22 Jan 90, 29 Jul 97

x=floor(x(:));
xmin=min(x); xmax=max(x); xrange=xmax-xmin;
ind=( x*ones(1,xrange+1) )==( ones(size(x))*(xmin:xmax) );
ind=ind(:,any(ind));
if nargin>1
   ind(:,setto0)=[];
end;
