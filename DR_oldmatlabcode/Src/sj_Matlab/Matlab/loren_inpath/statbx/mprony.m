function [b,a,mu] = mprony(y,t,p,bstart);
%MPRONY Fits a sum of exponential functions to data by a modified Prony
%	method.
%
%	B = MPRONY(Y,T,P) approximates Y(i) by a weighted sum of P
%	exponentials  sum_j=1^p A(j)*exp(B(j)*T(i)).  The T(i) are
%	assumed equally spaced.  B and A may be complex, although the
%	sum will be real if the data are.
%
%	Starting values for the B(i) are optional supplied by
%	MPRONY(Y,T,P,BSTART).  The A(i) and the fitted values are
%	optionally returned by [B,A,MU] = MPRONY(Y,T,P,BSTART).
%
%	This function is very effective on smaller data sets with T(i+1)-T(i)
%	not too small and P not too large.

%	References:
%	Kahn, M., Mackisack, M. S., Osborne, M. R., and Smyth, G. K. (1992).
%	   On the consistency of Prony's method and related algorithms.
%	   J. Comput. Graph. Statist., 1, 329-349.
%	Osborne, M. R., and Smyth, G. K. (1995). A modified Prony algorithm
%	   for fitting sums of exponential functions.  SIAM J. Sci. Statist.
%	   Comput., 16, 119-138.

%	Gordon Smyth, U of Queensland, gks@maths.uq.oz.au
%	17 Dec 90.  Latest revision 29 Apr 98.

y = y(:); t = t(:);     % impose column structure
[n cy]=size(y);
dt=t(2)-t(1);

% form Y
Y=zeros(n-p,p+1);
i=1:(n-p);
for j=1:p+1,
   Y(:,j)=y(i+j-1);
end;

% starting values
if nargin < 4,
   [x d]=eig(Y'*Y); [l jmin]=min(diag(d)); c=x(:,jmin);
else
   c=poly(exp(bstart*dt)); c=c(:);
   c=c(p+1:-1:1)/norm(c);
end;

% form B
% X=zeros(n,n-p);
X=sparse( [],[],[],n,n-p,(n-p)*(p+1) );
for j=1:n-p,
   X(j:j+p,j)=conj(c);
end;
MY=(X'*X)\Y;
v=MY*c;
V=zeros(n,p+1);
for j=1:p+1,
   V(j:n-p+j-1,j)=v;
end;
B=( Y'*MY-V'*V )./(n-p);

% initial eigenvalues
lold=inf;
[x d]=eig(B); [l jmin]=min(diag(d)); c=x(:,jmin);
lfirst=l;

tol=10^(-15+log(n));
iter=0;
while (abs((l-lold)/lfirst) > 1e-5) & (rcond(B) > tol) & (iter<40);
   iter=iter+1;
   if iter==40;
      disp('MProny.  Max iterations reached.');
   end;
%  form B
   for j=1:n-p,
      X(j:j+p,j)=conj(c);
   end;
   MY=(X'*X)\Y;
   v=MY*c;
   V=zeros(n,p+1);
   for j=1:p+1,
      V(j:n-p+j-1,j)=v;
   end;
   B=( Y'*MY-V'*V )./(n-p);
%  inverse iteration
   lold=l;
   [x d]=eig(B); [l jmin]=min(diag(d)); c=x(:,jmin);
end;

% extract rate constants
b=log(roots( c(p+1:-1:1) ))/dt;

if nargout > 1,
   A=exp(t*b.');
   a=A\y;
   mu=A*a;
end;
