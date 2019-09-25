function f = poisp(x,lambda);
%POISP	POISP(X,LAMBDA) is the cumulative distribution function at X
%	of the Poisson distribution with mean LAMBDA.

% Reference: Johnson & Kotz, Discrete Distributions, p 114

% GKS 22 Sep 95

if lambda <= 0
   error('Poisson mean must be positive.')
elseif x < 0
   f = 0;
else
   f = 1 - gamma( floor(x+1), lambda );
end
