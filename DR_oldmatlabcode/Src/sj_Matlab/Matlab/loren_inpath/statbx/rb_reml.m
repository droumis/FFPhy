function [beta,sebeta,gamma,segamma]=rb_reml(y,x,z)
%RB_REML	Restricted maximum likelihood estimation for mixed linear models.
%	[BETA,SEBETA,GAMMA,SEGAMMA] = RB_REML(Y,X,Z) fits the model
%	Y = X*BETA + Z*U + E where BETA is fixed and U is
%	random.
%
%	GAMMA holds the variance components.  The errors E and
%	random effects U are assumed to have covariance matrices
%	EYE*GAMMA(1) and EYE*GAMMA(2) respectively.
%
%	SEBETA and SEGAMMA give standard errors for BETA and GAMMA
%	respectively.

%	GKS  19 Feb 94

% transform to independent observations
[mx,nx]=size(x);
[u,d,v]=svd(x);
c=u(:,nx+1:mx);
[u,d,v]=svd(c'*z);
dy=(u'*(c'*y)).^2;
dx=[ones(mx-nx,1) sum(d')'.^2];

% fit gamma glm identity link to dy with dx as covariates
lp=dy;
iter=0;
dev=1e6; devold=dev+1;
while abs(devold-dev) > 1e-8 & iter < 20,
   iter=iter+1;
   mu=lp;
   devold=dev;
   dev=2*sum( (dy-mu)./mu - log(dy./mu) );
   w=1.0./mu.^2; dz=(dy-mu)+lp;
   gamma=( dx'*((w*ones(1,2)).*dx) )\( dx'*(w.*dz) );
   gamma=max(gamma,0);
   lp=dx*gamma;
end;
if iter==20, disp('REML:  no convergence'); end;
phi=2;
segamma=sqrt(phi.*diag(inv( dx'*((w*ones(1,2)).*dx) )));

[u,d,v]=svd(z);
x=u'*x;
y=u'*y;
v=[ones(mx,1) sum(d')'.^2]*gamma;
info=x'*(x./(v*ones(1,nx)));
beta=info\( x'*(y./v) );
sebeta=sqrt(diag(inv(info)));
