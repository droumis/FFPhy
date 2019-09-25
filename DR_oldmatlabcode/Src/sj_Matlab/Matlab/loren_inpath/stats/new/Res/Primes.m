% PRIMES: find the first k prime numbers by brute force.
%
%   Usage: p = primes(k)
%
%       k = number of primes to be returned.
%       ------------------------------------
%       p = column vector of the first k primes.
%

function p = primes(k)
  p = [1;2];
  candidate = 2;

  while (k > length(p))
    candidate = candidate + 1;
    r = rem(candidate,(2:candidate-1));
    if (~any(r==0))                         % if (all(r>0))
      p = [p; candidate];
    end;
  end;
  
  return;