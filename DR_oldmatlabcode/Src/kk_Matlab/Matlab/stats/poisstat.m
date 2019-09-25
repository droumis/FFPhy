function [m,v]= poisstat(lambda);
%POISSTAT Mean and variance for the Poisson distribution.
%   [M,V] = POISSTAT(LAMBDA) returns the mean and variance of
%   the Poisson distribution with parameter LAMBDA.
%
%   See also POISSCDF, POISSFIT, POISSINV, POISSPDF, POISSRND.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.9.2.4 $  $Date: 2004/07/28 04:39:00 $

if nargin <  1, 
    error('stats:poisstat:TooFewInputs','Requires one input argument.'); 
end

% Initialize mean and variance to zero.
if isa(lambda,'single')
   m  = zeros(size(lambda),'single');
   v = zeros(size(lambda),'single');
else
   m  = zeros(size(lambda));
   v = zeros(size(lambda));
end

% Lambda must be positive.
k = find(lambda <= 0);
if any(k)
    m(k) = NaN;
    v(k) = NaN;
end

k = find(lambda > 0);
if any(k)
    m(k) = lambda(k);
    v(k) = lambda(k);
end
