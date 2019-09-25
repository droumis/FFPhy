function [g,g1,p,dev]=lik(y,x,beta,z,z1)
%LIK	Called by LOGIST.  Calculates likelihood for ordinal logistic
%	regression model.

e=exp( [z x]*beta ); e1=exp( [z1 x]*beta );
g=e./(1+e); g1=e1./(1+e1);
g=max( y==max(y),g ); g1=min( y>min(y),g1 );
p=g-g1;
dev=-2*sum(log(p));
