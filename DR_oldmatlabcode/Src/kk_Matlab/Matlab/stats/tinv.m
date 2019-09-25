function x = tinv(p,v);
%TINV   Inverse of Student's T cumulative distribution function (cdf).
%   X=TINV(P,V) returns the inverse of Student's T cdf with V degrees 
%   of freedom, at the values in P.
%
%   The size of X is the common size of P and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also TCDF, TPDF, TRND, TSTAT, ICDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.2

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 2.11.2.7 $  $Date: 2005/12/12 23:34:13 $

if nargin < 2, 
    error('stats:tinv:TooFewInputs','Requires two input arguments.'); 
end

[errorcode p v] = distchck(2,p,v);

if errorcode > 0
    error('stats:tinv:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero, or NaN for invalid d.f.
if isa(p,'single') || isa(v,'single')
    x = NaN(size(p),'single');
else
    x = NaN(size(p));
end

% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
x(p==0 & v > 0) = -Inf;
x(p==1 & v > 0) = Inf;

k0 = (0<p & p<1) & (v > 0);

% Invert the Cauchy distribution explicitly
k = find(k0 & (v == 1));
if any(k)
  x(k) = tan(pi * (p(k) - 0.5));
end

% For small d.f., call betainv which uses Newton's method
k = find(k0 & (v < 1000));
if any(k)
    q = p(k) - .5;
    df = v(k);
    z = betainv(1-2.*abs(q),df/2,0.5);
    x(k) = sign(q) .* sqrt(df ./ z - df);
end

% For large d.f., use Abramowitz & Stegun formula 26.7.5
% k = find(p>0 & p<1 & ~isnan(x) & v >= 1000);
k = find(k0 & (v >= 1000));
if any(k)
   xn = norminv(p(k));
   df = v(k);
   x(k) = xn + (xn.^3+xn)./(4*df) + ...
           (5*xn.^5+16.*xn.^3+3*xn)./(96*df.^2) + ...
           (3*xn.^7+19*xn.^5+17*xn.^3-15*xn)./(384*df.^3) +...
           (79*xn.^9+776*xn.^7+1482*xn.^5-1920*xn.^3-945*xn)./(92160*df.^4);
end