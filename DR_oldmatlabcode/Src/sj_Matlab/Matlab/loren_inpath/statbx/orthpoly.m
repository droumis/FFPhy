function [op]=orthpoly(x,n,c)
%ORTHPOLY  ORTHPOLY(X,N) calculates the orthogonal polynomials up
%	to order N corresponding to vector X.
%
%	ORTHPOLY(X,N,C) calculates constrained polynomials, so that
%	the constant vector is not included.

%	GKS 5 Mar 92. Revised 5 Dec 92.

if n>=0;
   x=x(:);
   if nargin<3,
      op(:,1)=ones(size(x));
   else
      op(:,1)=x;
   end;
   if n>=1,
      if nargin<3,
         op(:,2)=x-mean(x);
      else
         op(:,2)=x.*(x-sum(x.^3)./sum(x.^2));
      end;
      for j=3:n+1-(nargin==3);
         a=sum(x.*op(:,j-1).^2)./sum(op(:,j-1).^2);
         b=sum(x.*op(:,j-1).*op(:,j-2))./sum(op(:,j-2).^2);
         op(:,j)=(x-a).*op(:,j-1)-b.*op(:,j-2);
      end;
   end;
else
   op=[];
end;
