function q = normq(p);
%NORMQ	NORMQ(P) is the inverse cumulative normal distribution function.

%	GKS 8 Oct 93

q = erfinv(2*p-1) * sqrt(2);
