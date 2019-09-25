function p = randchoose(n,m)
% randchoose - return a random choice.
%
%	RANDCHOOSE(N, M) returns a random choice of M numbers from 1:N.
%

if (m > n)
  error('m must be less than n');
end

p = [1:m] + floor(sort(rand(1,m))*(n-m+1));

