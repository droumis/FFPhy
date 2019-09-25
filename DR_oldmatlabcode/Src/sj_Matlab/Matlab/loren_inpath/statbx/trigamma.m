function t=trigamma(x);
%TRIGAMMA TRIGAMMA(X) is the trigamma function at X.

t=x;
for i=1:prod(size(x));
  t(i)=trigamms(x(i));
end;
