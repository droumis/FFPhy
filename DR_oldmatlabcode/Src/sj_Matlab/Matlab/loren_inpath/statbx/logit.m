function y=logit(x);
%LOGIT	LOGIT(X) = LOG(X./(1-X))

y=log(x./(1-x));
