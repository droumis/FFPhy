function [beta,sebeta,gamma,segamma,dev]=rb_ml(y,x,z)
%RB_ML	Maximum likelihood estimation for randomized block experiments.
%	[BETA,SEBETA,GAMMA,SEGAMMA,DEV] = RB_ML(Y,X,Z) fits the model
%	Y = X*BETA + Z*U + E where BETA is fixed and U is
%	random.
%
%	GAMMA holds the variance components.  The errors E and
%	random effects GAMMA are assumed to have covariances
%	EYE*GAMMA(1) and EYE*GAMMA(2) respectively.
%
%	SEBETA and SEGAMMA give standard errors.  DEV is minus twice
%	the log-likelihood.

%	GKS  19 Feb 94

% transform to independent observations
[mx,nx]=size(x);
[u,d,v]=svd(z);
y=u'*y;
x=u'*x;
z=[ones(mx,1) sum(d')'.^2];

% fit double generalized linear model
const=mx*log(2*pi);
iter=0;
mu=x*(x\y);
d=(y-mu).^2;
phi=z*(z\d);
devold=sum(log(phi))+sum(d./phi)+const;
dev=devold+1;
while abs(devold-dev) > 1e-6,

   % mean submodel: normal
   beta=( x./(sqrt(phi)*ones(1,nx)) )\( y./sqrt(phi) );
   mu=x*beta;
   d=(y-mu).^2;

   % dispersion submodel: gamma identity-link
   zd=d;
   gamma=( z./(phi*ones(1,2)) )\( zd./phi );
   gamma=max(gamma,0);
   phi=z*gamma;

   dev=sum(log(phi))+sum(d./phi)+const;
   devold=dev;
   iter=iter+1;
end;
if iter==20, disp('ML:  no convergence'); end;

sebeta=sqrt(diag(inv(x'*(x./(phi*ones(1,nx))))));
segamma=sqrt(2*diag(inv(z'*(z./(phi.^2*ones(1,2))))));
