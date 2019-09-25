function [dl,d2l]=derivs(x,z,z1,g,g1,p)
%DERIVS	Called by LOGIST.  Calculates derivates of log-likelihood for
%	ordinal logistic regression model.

%	Calls DMULT.

% first derivative
v=g.*(1-g)./p; v1=g1.*(1-g1)./p;
dlogp=[dmult(v,z)-dmult(v1,z1) dmult(v-v1,x)];
dl=sum(dlogp)';

% second derivative
w=v.*(1-2*g); w1=v1.*(1-2*g1);
d2l=[z x]'*dmult(w,[z x])-[z1 x]'*dmult(w1,[z1 x])-dlogp'*dlogp;
