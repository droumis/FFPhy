function y=gammad(x,alpha,beta);
%GAMMAD	y=GAMMAD(x,alpha,beta) is the density of a Gamma(alpha,beta)
%	random variable at x.

%	GKS 15 May 92

y=x.^(alpha-1).*exp(-x./beta)./gamma(alpha)./(beta.^alpha);
