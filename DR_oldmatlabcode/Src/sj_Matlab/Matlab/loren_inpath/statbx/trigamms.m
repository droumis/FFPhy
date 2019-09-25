function t=trigamms(x);
%TRIGAMMS TRIGAMMS(X) is the trigamma function for scalar argument X.
%	From Algorithm AS121 (Schneider, 1978).
%
%	GKS  Last revised 13 Jan 98
a=1e-8; b=12;
b2=1/6; b4=-1/30; b6=1/42;
if x<=0;
   disp('trigamma arg <= 0');
   t=0;
elseif x<a;
   t=1./(x.*x);
elseif x>b;
   y=1./(x.*x);
   t=.5*y+(1+y.*(b2+y.*(b4+y.*(b6+y.*b4))))./x;
else;
   t=1./(x.*x)+trigamms(x+1);
end;
