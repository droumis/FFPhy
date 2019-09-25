function d=digammas(x);
%DIGAMMAS DIGAMMAS(X) is the digamma function for scalar argument X.
%	From Algorithm AS103 (Bernardo, 1976).
%
%	GKS  Last revised 13 Jan 98
   a=1e-8; b=10;
   s3=1/12; s4=1/120; s5=1/252; g=-0.57721566490153;
   if x<=0,
      disp('digamma argument <= 0');
      d=0;
   elseif x<a,
      d=g-1./x;
   elseif x>b;
      y=1./(x.*x);
      d=log(x)-0.5./x-y.*(s3-y.*(s4-y.*s5));
   else;
      d=digammas(x+1)-1./x;
   end;
