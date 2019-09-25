% PERMNEXT: Given a permutation, finds the lexicographically next
%           permutation in the list of N! possible permutations.  
%           Wraps around from last to first permutation.
%
%       Usage: next = permnext(current)
%
%           current = permuted vector of integers 1:n
%           next =    next vector
%

% Based on Tucker,A. 1995. Applied Combinatorics, pp. 225-226. Wiley.

% RE Strauss, 2/25/95

function next = permnext(current)
  n = length(current);
  next = current;

  i = n-1;
  while (next(i) > next(i+1))        % Find rightmost switchable pair
    i = i-1;
    if (i<1)
      next(1:n) = next(n:-1:1);
      return;
    end;
  end;

  h = n;
  while (next(i) > next(h))           % Find rightmost h such that A(h) > A(i)
    h = h-1;
  end;
  temp = next(i);                     % Swap A(h),A(i)
  next(i) = next(h);
  next(h) = temp;

  next(i+1:n) = next(n:-1:i+1);       % Reverse rightmost remaining sequence
  return;

