function y = hygepdf(x,m,k,n)
%HYGEPDF Hypergeometric probability density function.
%   Y = HYGEPDF(X,M,K,N) returns the hypergeometric probability
%   density function at X with integer parameters M, K, and N.
%   Note: The density function is zero unless X is an integer.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also HYGECDF, HYGEINV, HYGERND, HYGESTAT, PDF.

%   Reference:
%      [1]  Mood, Alexander M., Graybill, Franklin A. and Boes, Duane C.,
%      "Introduction to the Theory of Statistics, Third Edition", McGraw Hill
%      1974 p. 91.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.11.2.5 $  $Date: 2005/11/18 14:28:05 $

if nargin < 4,
    error('stats:hygepdf:TooFewInputs','Requires four input arguments.');
end

[errorcode x m k n] = distchck(4,x,m,k,n);

if errorcode > 0
    error('stats:hygepdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
if isa(x,'single') || isa(m,'single') || isa(k,'single') || isa(n,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end

% Return NaN for values of the parameters outside their respective limits.
k1 = (m < 0 | k < 0 | n < 0 | round(m) ~= m | round(k) ~= k ...
                    | round(n) ~= n | n > m | k > m);
y(k1) = NaN;

% Remove values of X for which Y is zero by inspection.
k2 = (m - k - n + x + 1 <= 0 | x < 0 | round(x) ~= x | x > n | x > k);

% Find integer values of x that are within the correct range
kc = ~(k1 | k2);
if any(kc),
    if ~isfloat(x), x = double(x); end
    if ~isfloat(m), m = double(m); end
    if ~isfloat(k), k = double(k); end
    if ~isfloat(n), n = double(n); end
    kx = gammaln(k(kc)+1)-gammaln(x(kc)+1)-gammaln(k(kc)-x(kc)+1);
    mn = gammaln(m(kc)+1)-gammaln(n(kc)+1)-gammaln(m(kc)-n(kc)+1);
    mknx = gammaln(m(kc)-k(kc)+1)-gammaln(n(kc)-x(kc)+1)-gammaln(m(kc)-k(kc)-(n(kc)-x(kc))+1);
    y(kc) = exp(kx + mknx - mn);
end
