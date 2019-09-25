function d=digamma(x);
%DIGAMMA DIGAMMA(X) is the digamma function at X.

d=x;
for i=1:prod(size(x));
  d(i)=digammas(x(i));
end;
