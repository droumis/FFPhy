function p = chebyshev_poly(n)
%CHEBYSHEV_POLY Coefficients of the Chebyshev polynomial of the first kind.
%
%   P = CHEBYSHEV_POLY(N) returns the coefficients of the Nth-degree Chebyshev
%   polynomial of the first kind. P is a row vector of length N+1 containing
%   the polynomial coefficients in descending powers (suitable for evaluation
%   with POLYVAL). N must be a real finite non-negative integer.
%
%Written by SMK, 2010 January 30.
%

if ~isscalar(n) || ~isreal(n) || ~isfinite(n) || (n < 0) || (round(n) ~= n)
  error('N must be a real finite non-negative integer');
end

if (n == 0)
  p = 1;
else
  p = zeros([1, n+1]);
  q = p;
  p(end-1) = 1;
  q(end) = 1;
  for i = 2:n
    r = 2.*circshift(p,[1 -1]) - q;
    q = p;
    p = r;
  end
end

