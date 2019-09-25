function r=runmed(x,wind)
%RUNMED	RUNMED(X) is the running median of 3 of X.
%	If X is a matrix, its components are taken in column order.
%	RUNMED(X,WIND) is the running median with window width WIND.

%	GKS  6 January 93

if nargin<2, wind=3; end;
x=x(:); [m n]=size(x);

X=zeros(m-wind+1,wind);
for c=1:wind,
  X(:,c)=x(c:m-wind+c);
end;
r=median(X')';
