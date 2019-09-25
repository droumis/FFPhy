function [v,l] = invits(A,B,lstart,vstart)
%INVITS	[v,l]=INVITS(A) finds the eigenvalue l of A closest to
%       zero and its eigenvector v by inverse iteration.
%	This can be faster than EIG if the eigenvalue nearest zero
%	is well separated from the other eigenvalues.
%	INVITS(A,B) solves the generalized eigenproblem A*v=l*B*v
%	in which B is Hermitian but may be singular.
%	INVITS(A,B,l,v) provides starting values for l and v.

%	GKS  18 May 90.  Modified  13 Dec 90, 21 Feb 91, 2 Jun 93, 27 Apr 98.
%	Should rewrite to avoid rcond calls.

% check input
[ma,na]=size(A);
if nargin == 1
   B=eye(ma);
else
   if any(any(B ~= B'))
      error('INVITS Second matrix must be Hermitian.')
   end;
end;
[mb,nb]=size(B);
if std([ma na mb nb]) > 0
   error('INVITS Matrices must be square and of equal dimension.')
end;

% starting values
if nargin >= 3
   l=lstart;
else
   l=0;
end;
if nargin >= 4
   v=vstart;
else 
   [v,imin]=min(sum(A.^2));
   v=zeros(na,1);
   v(imin)=1;
end;

% ensure at least one iteration
if rcond(A) < 4e-15,
   l=norm(A,inf)/norm(B,inf)*1e-14;
end;

% inverse iteration
iter=0;
while rcond(A-l*B) > 5e-15
   iter=iter+1;
   if iter>50
      rcond(A-l*B)
      error('INVITS  Too many iterations.')
   end;
   w=(A-l*B)\(B*v);
   Bw=B*w;
   s=w'*Bw;
   dl=(v'*Bw)/s;
   l=l+dl;
   v=w/sqrt(s);
end;
