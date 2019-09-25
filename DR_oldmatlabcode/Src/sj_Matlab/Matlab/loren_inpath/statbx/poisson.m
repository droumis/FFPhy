function [beta,mu,dev,df,se]=poisson(y,x,offset,print);
%POISSON Fit Poisson generalized linear model with log-link.
%	[BETA,MU,DEV,DF,SE]=POISSON(Y,X,OFFSET,PRINT)
%	All input and output arguments except Y are optional.
%
%	Y - response vector
%	X - matrix of covariates, including the constant vector if required
%	OFFSET - offset if required
%	PRINT - enter argument if output required each iteration
%
%	BETA - regression parameter estimates
%	SE - associated standard errors
%	MU - fitted values
%	DEV - residual deviance
%	DF - residual degrees of freedom

% GKS 28 March 95

% Initialize
y=y(:);
[my ny]=size(y);
if nargin<2, x=ones(my,ny); end;
[mx nx]=size(x);
if nargin<3, offset=0; end;

% Starting values
yadj=y+0.5.*(y==0);
lp=log(y+0.5)-offset;

% Iteratively reweighted least squares
dev=1e6; devold=dev+1;
if nargin==4, disp('dev'); end;
while abs(devold-dev) > 1e-8;
   mu=exp(offset+lp);
   devold=dev;
   dev=2*sum( y.*log(yadj./mu) - (y-mu) );
   z=(y-mu)./mu+lp;
   beta=( x'*((mu*ones(1,nx)).*x) )\( x'*(mu.*z) );
   lp=x*beta;
   if nargin==4, disp(dev); end;
end;

if nargout>3, df=my-nx; end;
if nargout>4, se=sqrt(diag(inv( x'*((mu*ones(1,nx)).*x) ))); end;
