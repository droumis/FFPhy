function p = hygecdf(x,m,k,n)
%HYGECDF Hypergeometric cumulative distribution function.
%   P = HYGECDF(X,M,K,N) returns the hypergeometric cumulative
%   distribution function with parameters M, K, and N
%   at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also HYGEINV, HYGEPDF, HYGERND, HYGESTAT, CDF.

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.11.4.6 $  $Date: 2005/11/18 14:28:03 $

if nargin < 4,
    error('stats:hygecdf:TooFewInputs','Requires four input arguments.');
end

[errorcode x m k n] = distchck(4,x,m,k,n);

if errorcode > 0
    error('stats:hygecdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

%Initialize P to zero.
if isa(x,'single') || isa(m,'single') || isa(k,'single') || isa(n,'single')
   p = zeros(size(x),'single');
else
   p = zeros(size(x));
end


% Handle values of X for which P is zero by inspection.
k1 = (m - k - n + x + 1 <= 0 | x < 0);

% Handle values of X for which P is one by inspection.
k2 = (x >= n | x >= k);
p(k2) = 1;

% Return NaN for values of the parameters outside their respective limits.
k3 = (m < 0 | k < 0 | n < 0 | round(m) ~= m | round(k) ~= k | round(n) ~= n | ...
      n > m | k > m);
p(k3) = NaN;

kc = find(~(k1|k2|k3));

% Compute p when xx >= 0.
if any(kc)
    xx = floor(x);
    val = min([max(max(k(kc))) max(max(xx(kc))) max(max(n(kc)))]);
    i1 = (0:val)';
    compare = i1(:,ones(size(kc)));
    index = xx(kc);
	index = index(:);
    index = index(:,ones(size(i1)))';
    mbig = m(kc);
	mbig = mbig(:);
    mbig = mbig(:,ones(size(i1)))';
    kbig = k(kc);
	kbig = kbig(:);
    kbig = kbig(:,ones(size(i1)))';
    nbig = n(kc);
	nbig = nbig(:);
    nbig = nbig(:,ones(size(i1)))';
    p0 = hygepdf(compare,mbig,kbig,nbig);
    indicator = find(compare > index);
    p0(indicator) = 0;
    p(kc) = sum(p0,1);
end

% Make sure that round-off errors never make P greater than 1.
p(p > 1) = 1;
