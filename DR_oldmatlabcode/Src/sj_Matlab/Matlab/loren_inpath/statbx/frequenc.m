function [omega,rho,phi,mu] = frequenc(y,p,t);
%FREQUENC estimates a sum of sinusoidal signals.
%
%	The data Y is assumed to be a sum of P sinusoid signals
%	rho * sin(t*omega + phi) plus Gaussian noise.
%
%	[OMEGA,RHO,PHI,MU] = FREQUEN(Y,P,T)  returns the frequencies OMEGA,
%	amplitudes RHO, phases PHI and fitted values MU.
%
%	P defaults to one, and T defaults to 1:length(Y).

%	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
%	28 September 98.

y = y(:);    % impose column structure
[n cy]=size(y);
if nargin<2, p=1; end;
if nargin<3, t=(1:n)'; end;
p=2*p;       % Number of complex exponentials
X=sparse( [],[],[],n,n-p,(n-p)*(p+1) );

% form Q
Q=eye(p+1);
Q=Q+Q(p+1:-1:1,:);
Q(:,floor(p/2+2):p+1)=[];
Q=sqrt(Q/2);

% form Y
Y=zeros(n-p,p+1);
for j=1:p+1
   Y(:,j)=y(j:n-p+j-1);
end;
YQ=Y*Q;

% Constrained Pisarenko
B=( YQ'*YQ )./(n-p);
[x d]=eig(B); [l jmin]=min(diag(d)); g=x(:,jmin);
c=Q*g;

% Constrained ORA
gold=zeros(p/2+1,1);
iter=0;
while norm(g-gold) > 1e-6 & iter < 50;
   iter=iter+1;
   gold=g;
   for j=1:n-p,
      X(j:j+p,j)=conj(c);
   end;
   B=( YQ'*((X'*X)\YQ) )./(n-p);
   [x d]=eig(B); [l jmin]=min(diag(d)); g=x(:,jmin);
   c=Q*g;
end;

% Constrained least squares
gold=zeros(p/2+1,1);
iter=0;
while norm(g-gold) > 1e-6 & iter < 50;
   iter=iter+1;
   gold=g;
   for j=1:n-p,
      X(j:j+p,j)=conj(c);
   end;
   MY=(X'*X)\Y;
   v=MY*c;
   V=zeros(n,p+1);
   for j=1:p+1,
      V(j:n-p+j-1,j)=v;
   end;
   B=Q'*( Y'*MY-V'*V )*Q./(n-p);
   [x d]=eig(B); [l jmin]=min(diag(d)); g=x(:,jmin);
   c=Q*g;
end;

% Extract frequencies
om=sort(imag(log(roots( c(p+1:-1:1) ))));
omega=om(p/2+1:p);

% Extract amplitudes and phases
if nargout>1,
   A=[sin(t*omega.') cos(t*omega.')];
   b=A\y;
   b1=b(1:p/2);
   b2=b(p/2+1:p);
   phi=atan(b2./b1);
   rho=sqrt(b1.*b1+b2.*b2);
end;

% Extract fitted values
if nargout>3, mu=A*b; end;
