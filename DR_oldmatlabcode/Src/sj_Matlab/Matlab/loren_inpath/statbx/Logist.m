function [beta,mu,dev,df,se]=logist(y,n,x,offset,print);
%LOGIST Fit logistic regression model.
%	[BETA,MU,DEV,DF,SE]=LOGIST(Y,N,X,OFFSET,PRINT)
%	All input and output arguments except Y are optional.
%
%	Y - response vector containing binomial counts
%	N - number of trials for each count. Y is assumed to be binomial(p,N).
%	X - matrix of covariates, including the constant vector if required
%	OFFSET - offset if required
%	PRINT - enter any argument if output required each iteration
%
%	BETA - regression parameter estimates
%	SE - associated standard errors
%	MU - fitted values
%	DEV - residual deviance
%	DF - residual degrees of freedom

% GKS 18 May 2002

% Initialize
y=y(:);
[my ny]=size(y);
if nargin<2, disp('Must specify binomial N'); return; end;
n=n(:);
if (min(n-y)<0), disp('Binomial counts Y must be less than or equal to binomial N'); return; end;
if nargin<3, x=ones(my,ny); end;
[mx nx]=size(x);
if nargin<4, offset=0; end;

% Starting values
y0=y+0.5.*(y==0);
yn=y-0.5.*(y==n);
lp=log((y+0.5)./(n-y+0.5))-offset;

% Iteratively reweighted least squares
dev=1e6; devold=dev+1;
if nargin==5, disp('dev'); end;
while abs(devold-dev) > 1e-8;
   p=exp(offset+lp);
   p=p./(1+p);
   mu=n.*p;
   devold=dev;
   dev=2.*sum(y.*log(y0./mu) + (n-y).*log((n-yn)./(n-mu)));
   v=mu.*(n-mu)./n;
   z=(y-mu)./v+lp;
   beta=( x'*((v*ones(1,nx)).*x) )\( x'*(v.*z) );
   lp=x*beta;
   if nargin==5, disp(dev); end;
end;

if nargout>3, df=my-nx; end;
if nargout>4, se=sqrt(diag(inv( x'*((v*ones(1,nx)).*x) ))); end;
